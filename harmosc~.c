/*
	harmosc~ - Computationally efficient harmonic oscillator for Max/MSP
	
	Features:
	- Variable harmonic count at instantiation
	- Dynamic falloff control for harmonic amplitude distribution  
	- Selective harmonic activation (all/odd/even)
	- Wavetable-based synthesis for performance
	
	Usage: [harmosc~ <freq> <harmonics> <falloff> <detune>]
	Defaults: 440 Hz, 8 harmonics, 0.0 falloff, 0.0 detune
	All arguments optional
*/

#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"
#include <math.h>

#define TABLE_SIZE 4096
#define PI 3.14159265358979323846

typedef struct _harmosc {
    t_pxobject x_obj;           // MSP object header
    
    // Oscillator state
    double phase;               // Master phase accumulator (0-1)
    double freq;                // Fundamental frequency
    double sr;                  // Sample rate
    double sr_recip;            // 1/sr for efficiency
    
    // Harmonic control
    int num_harmonics;          // Total number of harmonics
    double falloff;             // Falloff parameter (-1 to 1)
    double detune;              // Detune parameter (0-1)
    double *amplitudes;         // Pre-calculated amplitude array
    char *harmonic_states;      // On/off state for each harmonic
    double *detune_offsets;     // Random detune offsets for each harmonic
    
    // Optimization
    double *phase_increments;   // Pre-calculated phase increments
    double *sine_table;         // Wavetable for efficiency
    int table_size;             // Size of sine table
    int table_mask;             // For fast modulo operations
    
    // State flags
    char all_on;                // All harmonics active
    char odd_only;              // Only odd harmonics
    char even_only;             // Only even harmonics
    char custom_amps;           // Using custom amplitude values
} t_harmosc;

// Function prototypes
void *harmosc_new(t_symbol *s, long argc, t_atom *argv);
void harmosc_free(t_harmosc *x);
void harmosc_assist(t_harmosc *x, void *b, long m, long a, char *s);
void harmosc_float(t_harmosc *x, double f);
void harmosc_dsp64(t_harmosc *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void harmosc_perform64(t_harmosc *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);

// Message handlers
void harmosc_falloff(t_harmosc *x, double f);
void harmosc_detune(t_harmosc *x, double d);
void harmosc_all(t_harmosc *x);
void harmosc_odd(t_harmosc *x);
void harmosc_even(t_harmosc *x);
void harmosc_amps(t_harmosc *x, t_symbol *s, long argc, t_atom *argv);

// Internal functions
void harmosc_build_sine_table(t_harmosc *x);
void harmosc_update_phase_increments(t_harmosc *x);
void harmosc_calculate_amplitudes(t_harmosc *x);
void harmosc_generate_detune_offsets(t_harmosc *x);

// Global class pointer
static t_class *harmosc_class;

void ext_main(void *r) {
    t_class *c;
    
    c = class_new("harmosc~", (method)harmosc_new, (method)harmosc_free, 
                  sizeof(t_harmosc), NULL, A_GIMME, 0);
    
    class_addmethod(c, (method)harmosc_dsp64, "dsp64", A_CANT, 0);
    class_addmethod(c, (method)harmosc_assist, "assist", A_CANT, 0);
    class_addmethod(c, (method)harmosc_float, "float", A_FLOAT, 0);
    
    // Message handlers
    class_addmethod(c, (method)harmosc_falloff, "falloff", A_FLOAT, 0);
    class_addmethod(c, (method)harmosc_detune, "detune", A_FLOAT, 0);
    class_addmethod(c, (method)harmosc_all, "all", 0);
    class_addmethod(c, (method)harmosc_odd, "odd", 0);
    class_addmethod(c, (method)harmosc_even, "even", 0);
    class_addmethod(c, (method)harmosc_amps, "amps", A_GIMME, 0);
    
    class_dspinit(c);
    class_register(CLASS_BOX, c);
    harmosc_class = c;
}

void *harmosc_new(t_symbol *s, long argc, t_atom *argv) {
    t_harmosc *x = (t_harmosc *)object_alloc(harmosc_class);
    
    if (x) {
        dsp_setup((t_pxobject *)x, 0);  // No signal inlets
        outlet_new(x, "signal");         // One signal outlet
        
        // Parse arguments: [fundamental freq] [n harmonics] [falloff] [detune]
        x->num_harmonics = 8;   // Default
        x->freq = 440.0;        // Default  
        x->falloff = 0.0;       // Default
        x->detune = 0.0;        // Default
        
        // Argument 1: Fundamental frequency
        if (argc >= 1 && (atom_gettype(argv) == A_FLOAT || atom_gettype(argv) == A_LONG)) {
            double freq_arg = atom_gettype(argv) == A_FLOAT ? atom_getfloat(argv) : (double)atom_getlong(argv);
            x->freq = CLAMP(freq_arg, 0.1, 20000.0);
        }
        
        // Argument 2: Number of harmonics
        if (argc >= 2 && atom_gettype(argv + 1) == A_LONG) {
            x->num_harmonics = CLAMP(atom_getlong(argv + 1), 1, 64);
        }
        
        // Argument 3: Falloff
        if (argc >= 3 && (atom_gettype(argv + 2) == A_FLOAT || atom_gettype(argv + 2) == A_LONG)) {
            double falloff_arg = atom_gettype(argv + 2) == A_FLOAT ? atom_getfloat(argv + 2) : (double)atom_getlong(argv + 2);
            x->falloff = CLAMP(falloff_arg, -1.0, 1.0);
        }
        
        // Argument 4: Detune
        if (argc >= 4 && (atom_gettype(argv + 3) == A_FLOAT || atom_gettype(argv + 3) == A_LONG)) {
            double detune_arg = atom_gettype(argv + 3) == A_FLOAT ? atom_getfloat(argv + 3) : (double)atom_getlong(argv + 3);
            x->detune = CLAMP(detune_arg, 0.0, 1.0);
        }
        
        // Initialize oscillator state
        x->phase = 0.0;
        x->sr = sys_getsr();
        x->sr_recip = 1.0 / x->sr;
        
        // Initialize state flags
        x->all_on = 1;
        x->odd_only = 0;
        x->even_only = 0;
        x->custom_amps = 0;
        
        // Allocate arrays
        x->amplitudes = (double *)sysmem_newptr(x->num_harmonics * sizeof(double));
        x->harmonic_states = (char *)sysmem_newptr(x->num_harmonics * sizeof(char));
        x->phase_increments = (double *)sysmem_newptr(x->num_harmonics * sizeof(double));
        x->detune_offsets = (double *)sysmem_newptr(x->num_harmonics * sizeof(double));
        
        if (!x->amplitudes || !x->harmonic_states || !x->phase_increments || !x->detune_offsets) {
            object_error((t_object *)x, "Failed to allocate memory");
            return NULL;
        }
        
        // Initialize harmonic states (all on)
        for (int i = 0; i < x->num_harmonics; i++) {
            x->harmonic_states[i] = 1;
        }
        
        // Build sine table
        harmosc_build_sine_table(x);
        
        // Generate detune offsets and calculate initial values
        harmosc_generate_detune_offsets(x);
        harmosc_update_phase_increments(x);
        harmosc_calculate_amplitudes(x);
        
        post("harmosc~: Created with %.1f Hz, %d harmonics, falloff %.2f, detune %.2f", 
             x->freq, x->num_harmonics, x->falloff, x->detune);
    }
    
    return x;
}

void harmosc_free(t_harmosc *x) {
    if (x->sine_table) {
        sysmem_freeptr(x->sine_table);
    }
    if (x->amplitudes) {
        sysmem_freeptr(x->amplitudes);
    }
    if (x->harmonic_states) {
        sysmem_freeptr(x->harmonic_states);
    }
    if (x->phase_increments) {
        sysmem_freeptr(x->phase_increments);
    }
    if (x->detune_offsets) {
        sysmem_freeptr(x->detune_offsets);
    }
    dsp_free((t_pxobject *)x);
}

void harmosc_assist(t_harmosc *x, void *b, long m, long a, char *s) {
    if (m == ASSIST_INLET) {
        sprintf(s, "(float) Frequency in Hz\n(amps) Custom amplitude list\n(falloff) Falloff parameter\n(detune) Detune amount\n(all/odd/even) Harmonic selection");
    } else {
        sprintf(s, "(signal) Harmonic oscillator output");
    }
}

void harmosc_float(t_harmosc *x, double f) {
    // Frequency input
    x->freq = CLAMP(f, 0.1, 20000.0);
    harmosc_update_phase_increments(x);
}

void harmosc_build_sine_table(t_harmosc *x) {
    x->table_size = TABLE_SIZE;
    x->table_mask = TABLE_SIZE - 1;
    x->sine_table = (double *)sysmem_newptr(TABLE_SIZE * sizeof(double));
    
    if (!x->sine_table) {
        object_error((t_object *)x, "Failed to allocate sine table");
        return;
    }
    
    for (int i = 0; i < TABLE_SIZE; i++) {
        x->sine_table[i] = sin(2.0 * PI * i / TABLE_SIZE);
    }
}

void harmosc_update_phase_increments(t_harmosc *x) {
    double base_increment = x->freq * x->sr_recip;
    
    for (int i = 0; i < x->num_harmonics; i++) {
        double harmonic_frequency = i + 1;  // Base harmonic number (1, 2, 3, ...)
        
        // Apply detune: 0.0 = harmonic, 1.0 = ±50 cents maximum
        if (x->detune > 0.0) {
            // Convert cents to frequency ratio: 2^(cents/1200)
            double cents_offset = x->detune_offsets[i] * x->detune;
            double frequency_ratio = pow(2.0, cents_offset / 1200.0);
            double detuned_frequency = harmonic_frequency * frequency_ratio;
            x->phase_increments[i] = base_increment * detuned_frequency;
        } else {
            x->phase_increments[i] = base_increment * harmonic_frequency;
        }
    }
}

void harmosc_calculate_amplitudes(t_harmosc *x) {
    double total = 0.0;
    
    // Skip automatic calculation if using custom amplitudes
    if (x->custom_amps) {
        // Just apply harmonic states and normalize
        for (int i = 0; i < x->num_harmonics; i++) {
            x->amplitudes[i] *= x->harmonic_states[i];
            total += x->amplitudes[i];
        }
        
        // Normalize to maintain consistent output level
        if (total > 0.0) {
            for (int i = 0; i < x->num_harmonics; i++) {
                x->amplitudes[i] /= total;
            }
        }
        return;
    }
    
    // Calculate amplitudes with bipolar falloff behavior
    // -1.0 = only fundamental, 0.0 = equal harmonics, +1.0 = only highest harmonic
    for (int i = 0; i < x->num_harmonics; i++) {
        double harmonic_number = i + 1;  // 1-indexed harmonic number
        double highest_harmonic = x->num_harmonics;  // Highest harmonic number
        
        if (x->falloff == -1.0) {
            // Only fundamental
            x->amplitudes[i] = (i == 0) ? 1.0 : 0.0;
        } else if (x->falloff == 1.0) {
            // Only highest harmonic
            x->amplitudes[i] = (i == x->num_harmonics - 1) ? 1.0 : 0.0;
        } else if (x->falloff == 0.0) {
            // Equal amplitudes
            x->amplitudes[i] = 1.0;
        } else if (x->falloff < 0.0) {
            // Negative range: -1 to 0 (fundamental to equal)
            // Exponential decay from fundamental, getting weaker as we approach 0
            double decay_strength = -x->falloff;  // 0.0 to 1.0
            double decay_exponent = decay_strength * 3.0;  // Scale decay rate
            x->amplitudes[i] = pow(harmonic_number, -decay_exponent);
        } else {
            // Positive range: 0 to 1 (equal to highest harmonic)  
            // Exponential decay from highest harmonic, getting stronger as we approach 1
            double decay_strength = x->falloff;  // 0.0 to 1.0
            double decay_exponent = decay_strength * 3.0;  // Scale decay rate
            // Reverse the harmonic order for decay calculation
            double reverse_harmonic = highest_harmonic - harmonic_number + 1;
            x->amplitudes[i] = pow(reverse_harmonic, -decay_exponent);
        }
        
        // Apply harmonic states during calculation
        x->amplitudes[i] *= x->harmonic_states[i];
        total += x->amplitudes[i];
    }
    
    // Normalize to maintain consistent output level
    if (total > 0.0) {
        for (int i = 0; i < x->num_harmonics; i++) {
            x->amplitudes[i] /= total;
        }
    }
}

void harmosc_dsp64(t_harmosc *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags) {
    // Update sample rate if changed
    if (x->sr != samplerate) {
        x->sr = samplerate;
        x->sr_recip = 1.0 / samplerate;
        harmosc_update_phase_increments(x);
    }
    
    object_method(dsp64, gensym("dsp_add64"), x, harmosc_perform64, 0, NULL);
}

void harmosc_perform64(t_harmosc *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam) {
    double *out = outs[0];
    double phase = x->phase;
    double *sine_table = x->sine_table;
    int table_mask = x->table_mask;
    double base_increment = x->phase_increments[0];
    
    // Process audio block
    for (long i = 0; i < sampleframes; i++) {
        double sample = 0.0;
        
        // Additive synthesis loop - uses individual phase increments for detuning
        for (int h = 0; h < x->num_harmonics; h++) {
            if (x->amplitudes[h] > 0.0) {
                // Use individual phase increment ratio for detuning
                double phase_ratio = x->phase_increments[h] / base_increment;
                double harmonic_phase = phase * phase_ratio;
                // Wrap phase to 0-1 range
                harmonic_phase = harmonic_phase - floor(harmonic_phase);
                int table_index = (int)(harmonic_phase * TABLE_SIZE) & table_mask;
                sample += sine_table[table_index] * x->amplitudes[h];
            }
        }
        
        *out++ = sample;
        
        // Update phase
        phase += base_increment;
        if (phase >= 1.0) phase -= 1.0;
    }
    
    x->phase = phase;
}

// Message handlers
void harmosc_falloff(t_harmosc *x, double f) {
    x->falloff = CLAMP(f, -1.0, 1.0);
    x->custom_amps = 0;  // Clear custom amplitudes when using falloff
    harmosc_calculate_amplitudes(x);
    post("harmosc~: falloff set to %.3f", x->falloff);
}

void harmosc_detune(t_harmosc *x, double d) {
    x->detune = CLAMP(d, 0.0, 1.0);
    harmosc_update_phase_increments(x);
    post("harmosc~: detune set to %.3f", x->detune);
    
    // Debug: show actual cents detuning for first few harmonics
    for (int i = 1; i < MIN(4, x->num_harmonics); i++) {
        double cents = x->detune_offsets[i] * x->detune;
        post("  harmonic %d: %.1f cents (offset %.1f * detune %.3f)", 
             i+1, cents, x->detune_offsets[i], x->detune);
    }
}

void harmosc_all(t_harmosc *x) {
    x->all_on = 1;
    x->odd_only = 0;
    x->even_only = 0;
    
    for (int i = 0; i < x->num_harmonics; i++) {
        x->harmonic_states[i] = 1;
    }
    harmosc_calculate_amplitudes(x);
    post("harmosc~: all harmonics enabled");
}

void harmosc_odd(t_harmosc *x) {
    x->all_on = 0;
    x->odd_only = 1;
    x->even_only = 0;
    
    for (int i = 0; i < x->num_harmonics; i++) {
        // Always keep fundamental (i=0), then odd harmonics only
        if (i == 0) {
            x->harmonic_states[i] = 1;  // Fundamental always on
        } else {
            // Harmonic numbers are 1-indexed (fundamental = 1)
            x->harmonic_states[i] = ((i + 1) % 2 == 1) ? 1 : 0;
        }
    }
    harmosc_calculate_amplitudes(x);
    post("harmosc~: odd harmonics enabled (with fundamental)");
}

void harmosc_even(t_harmosc *x) {
    x->all_on = 0;
    x->odd_only = 0;
    x->even_only = 1;
    
    for (int i = 0; i < x->num_harmonics; i++) {
        // Always keep fundamental (i=0), then even harmonics only
        if (i == 0) {
            x->harmonic_states[i] = 1;  // Fundamental always on
        } else {
            // Harmonic numbers are 1-indexed (fundamental = 1)
            x->harmonic_states[i] = ((i + 1) % 2 == 0) ? 1 : 0;
        }
    }
    harmosc_calculate_amplitudes(x);
    post("harmosc~: even harmonics enabled (with fundamental)");
}

void harmosc_generate_detune_offsets(t_harmosc *x) {
    // Generate random detune offsets for each harmonic
    // Range: ±50 cents maximum (detune=1.0) for musical detuning
    for (int i = 0; i < x->num_harmonics; i++) {
        if (i == 0) {
            // Keep fundamental undetuned
            x->detune_offsets[i] = 0.0;
        } else {
            // Random offset between -50 and +50 cents
            double random_val = (double)rand() / RAND_MAX;  // 0.0 to 1.0
            x->detune_offsets[i] = (random_val * 100.0) - 50.0;  // -50 to +50 cents
            
            // Debug output to verify detuning amounts
            post("harmosc~: harmonic %d detune offset: %.1f cents", i+1, x->detune_offsets[i]);
        }
    }
}

void harmosc_amps(t_harmosc *x, t_symbol *s, long argc, t_atom *argv) {
    if (argc == 0) {
        object_error((t_object *)x, "amps: requires at least one amplitude value");
        return;
    }
    
    // Set custom amplitudes flag
    x->custom_amps = 1;
    
    // Clear state flags since we're using custom amplitudes
    x->all_on = 0;
    x->odd_only = 0;
    x->even_only = 0;
    
    // Process each amplitude value
    int num_values = MIN(argc, x->num_harmonics);
    for (int i = 0; i < num_values; i++) {
        double amp_value = 0.0;
        
        // Extract amplitude value from atom
        if (atom_gettype(argv + i) == A_FLOAT) {
            amp_value = atom_getfloat(argv + i);
        } else if (atom_gettype(argv + i) == A_LONG) {
            amp_value = (double)atom_getlong(argv + i);
        } else {
            object_error((t_object *)x, "amps: argument %d is not a number", i + 1);
            continue;
        }
        
        // Clamp amplitude to 0-1 range
        amp_value = CLAMP(amp_value, 0.0, 1.0);
        x->amplitudes[i] = amp_value;
        
        // Set harmonic state based on amplitude (0 = off, >0 = on)
        x->harmonic_states[i] = (amp_value > 0.0) ? 1 : 0;
    }
    
    // Set remaining harmonics to 0 if fewer values provided than harmonics
    for (int i = num_values; i < x->num_harmonics; i++) {
        x->amplitudes[i] = 0.0;
        x->harmonic_states[i] = 0;
    }
    
    // Apply normalization via calculate_amplitudes (but skip the automatic calculation)
    harmosc_calculate_amplitudes(x);
    
    post("harmosc~: custom amplitudes set for %d harmonics", num_values);
}