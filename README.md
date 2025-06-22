# harmosc~ - High-Performance Harmonic Oscillator for Max/MSP

A computationally efficient harmonic oscillator external that provides precise control over harmonic content through multiple methods: automatic falloff curves, harmonic filtering, and custom amplitude arrays.

## Features

- **Variable Harmonic Count**: Configure 1-64 harmonics at instantiation
- **Custom Amplitude Control**: Set individual harmonic amplitudes via list messages
- **Dynamic Falloff**: Bipolar falloff parameter for automatic harmonic distribution
- **Harmonic Filtering**: Selective activation of all/odd/even harmonics
- **Musical Detuning**: Randomized harmonic detuning for organic textures
- **Wavetable Synthesis**: Optimized sine table lookup for real-time performance
- **Multi-Argument Setup**: Flexible instantiation with optional parameters
- **Universal Binary**: Compatible with Intel and Apple Silicon Macs

## Usage

### Basic Instantiation
```
[harmosc~]                           // Default: 440 Hz, 8 harmonics
[harmosc~ 220]                       // Set frequency: 220 Hz, 8 harmonics
[harmosc~ 110 16]                    // Set frequency + harmonics: 110 Hz, 16 harmonics
[harmosc~ 440 12 -0.5]              // Add falloff: 440 Hz, 12 harmonics, fundamental emphasis
[harmosc~ 220 8 0.3 0.2]            // Full setup: 220 Hz, 8 harmonics, high emphasis, 20% detune
```

### Runtime Control
```
// Frequency control
[440(               // Set frequency to 440 Hz

// Falloff control (-1 to +1)
[falloff -1(        // Only fundamental
[falloff 0(         // Equal amplitudes
[falloff 1(         // Only highest harmonic

// Harmonic selection
[all(               // Enable all harmonics
[odd(               // Only odd harmonics (+ fundamental)
[even(              // Only even harmonics (+ fundamental)

// Detuning (0-1, adds ±50 cents maximum)
[detune 0.3(        // 30% of maximum detuning

// Custom amplitude arrays
[amps 1.0 0.5 0.3 0.8 0.2(          // Set specific harmonic levels
[amps 1.0 0.0 0.33 0.0 0.2(         // Square wave approximation
```

## Parameters

### Instantiation Arguments
All arguments are optional and can be provided in order:

1. **Frequency** (float, 0.1-20000 Hz, default: 440)
   - Fundamental frequency in Hz
   - Example: `[harmosc~ 220]`

2. **Harmonics** (int, 1-64, default: 8)
   - Number of harmonics to generate
   - Example: `[harmosc~ 440 16]`

3. **Falloff** (float, -1.0 to 1.0, default: 0.0)
   - Automatic amplitude distribution
   - -1.0 = only fundamental, 0.0 = equal, +1.0 = only highest
   - Example: `[harmosc~ 440 8 -0.5]`

4. **Detune** (float, 0.0-1.0, default: 0.0)
   - Random harmonic detuning amount
   - 0.0 = perfect harmonics, 1.0 = up to ±50 cents per harmonic
   - Example: `[harmosc~ 440 8 0.0 0.3]`

### Message Control

#### Falloff Parameter (-1.0 to 1.0)
Controls automatic harmonic amplitude distribution:
- **-1.0**: Only fundamental (pure sine wave)
- **-0.5**: Strong fundamental emphasis (exponential decay)
- **0.0**: Equal amplitudes (sawtooth-like spectrum)
- **+0.5**: High harmonic emphasis (reverse exponential)
- **+1.0**: Only highest harmonic

#### Custom Amplitude Arrays
```
[amps 1.0 0.5 0.3 0.8]    // Set harmonics 1-4 to specific levels
[amps 1.0 0.0 0.33]       // Harmonics 1,3 only (others = 0)
```
- Accepts 1 to N values (where N = number of harmonics)
- Values automatically clamped to 0.0-1.0 range
- Amplitudes are normalized to maintain consistent output level
- Overrides falloff parameter when used
- Remaining harmonics set to 0 if fewer values provided

#### Harmonic Selection
- **[all]**: Enable all harmonics (default state)
- **[odd]**: Enable odd harmonics only (1st, 3rd, 5th, etc.)
- **[even]**: Enable even harmonics only (2nd, 4th, 6th, etc.)
- **Note**: Fundamental (1st harmonic) always remains active

#### Detuning (0.0-1.0)
- **0.0**: Perfect harmonic ratios (1:2:3:4...)
- **0.5**: Moderate detuning (±25 cents average)
- **1.0**: Maximum detuning (±50 cents maximum)
- Each harmonic gets random offset within range
- Fundamental always remains undetuned

## Musical Applications

### Synthesis Techniques

**Additive Synthesis**
```
[harmosc~ 110 16]          // Rich harmonic content
[falloff -0.3(             // Natural decay curve
[detune 0.1(               // Slight organic detuning
```

**Formant-like Resonances**
```
[harmosc~ 440 12]
[amps 0.8 0.2 1.0 0.4 0.9 0.3 0.2 0.1]  // Emphasize specific harmonics
```

**Square Wave Synthesis**
```
[harmosc~ 220 16]
[odd(                      // Only odd harmonics
[falloff -0.2(             // Natural odd harmonic decay
```

**Sawtooth Wave Synthesis**
```
[harmosc~ 220 16]
[all(                      // All harmonics
[falloff 0.0(              // Equal amplitudes
```

**Spectral Animation**
```
[phasor~ 0.1]              // Slow LFO
[scale -1. 1. -1. 1.]      // Bipolar scaling
[falloff $1(               // Animate falloff parameter
```

### Advanced Techniques

**Dynamic Harmonic Content**
```
[metro 500]                // Regular updates
[drunk 8]                  // Random walk
[scale 0 16 0. 1.]         // Scale to amplitude range
[pak 1. 0. 0. 0. 0. 0. 0. 0.]  // Build amplitude list
[amps $1 $2 $3 $4 $5 $6 $7 $8( // Send to harmosc~
```

**Chord Generation**
```
[harmosc~ 220 8 -0.3]      // Root
[harmosc~ 275 8 -0.3]      // Fifth (5:4 ratio)
[harmosc~ 330 8 -0.3]      // Major third (3:2 ratio)
[*~ 0.7]                   // Mix levels
[+~]                       // Sum for chord
```

**Organic Textures**
```
[harmosc~ 55 16 0.2 0.4]   // Deep frequency, emphasis high, strong detune
[lores~ 2000 0.8]          // Filter for warmth
[*~ 0.8]                   // Control level
```

## Technical Implementation

### Performance Optimizations
- **Wavetable Lookup**: 4096-point sine table for efficient computation
- **Pre-calculated Values**: Phase increments and amplitudes cached
- **Conditional Processing**: Harmonics with zero amplitude skipped
- **Memory Efficiency**: Dynamic allocation based on harmonic count

### Algorithm Details
```c
// Core synthesis loop (simplified)
for (each_sample) {
    sample = 0.0;
    for (each_harmonic) {
        if (amplitude[h] > 0.0) {
            phase = master_phase * harmonic_ratio[h];
            sample += sine_table[phase] * amplitude[h];
        }
    }
    output = sample;
    master_phase += frequency_increment;
}
```

### Memory Usage
- **Base Object**: ~200 bytes
- **Per Harmonic**: ~40 bytes (amplitudes, states, phase increments, detune offsets)
- **Sine Table**: 32KB (shared, allocated once)
- **Total (16 harmonics)**: ~32.8KB per instance

## Build Instructions

### Requirements
- CMake 3.19 or later
- Xcode command line tools (macOS)
- Max SDK base (included)

### Building
```bash
cd source/audio/harmosc~
mkdir build && cd build
cmake -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
cmake --build .

# For M1/M2 Macs
codesign --force --deep -s - ../../../externals/harmosc~.mxo
```

### Verification
```bash
# Check universal binary
file ../../../externals/harmosc~.mxo/Contents/MacOS/harmosc~
# Should output: Mach-O universal binary with 2 architectures

# Test in Max
# Load harmosc~.maxhelp to verify functionality
```

## Files

- `harmosc~.c` - Main external implementation with A_GIMME list handling
- `CMakeLists.txt` - Build configuration
- `README.md` - This documentation
- `CLAUDE.md` - Development notes and implementation details
- `harmosc~.maxhelp` - Interactive help file

## Examples and Presets

### Waveform Approximations
```
// Square wave (odd harmonics)
[harmosc~ 220 16]
[amps 1.0 0.0 0.33 0.0 0.2 0.0 0.14 0.0 0.11 0.0 0.09]

// Sawtooth wave (all harmonics, 1/n decay)
[harmosc~ 220 16]  
[amps 1.0 0.5 0.33 0.25 0.2 0.16 0.14 0.125 0.11 0.1]

// Triangle wave (odd harmonics, 1/n² decay)
[harmosc~ 220 16]
[amps 1.0 0.0 0.11 0.0 0.04 0.0 0.02 0.0 0.012 0.0 0.008]
```

### Spectral Effects
```
// Comb filtering effect
[harmosc~ 110 16]
[amps 1.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 1.0]  // Every 4th harmonic

// Formant emphasis
[harmosc~ 220 16]
[amps 0.5 0.3 0.8 1.0 0.9 0.4 0.2 0.1]  // Emphasize 3rd-5th harmonics
```

## Compatibility

- **Max/MSP**: Version 8.0 or later
- **macOS**: 10.14 or later (universal binary)
- **Windows**: Not currently supported

## Development Patterns Demonstrated

This external showcases several Max SDK patterns:

1. **Multi-Argument Parsing**: Flexible instantiation with optional parameters
2. **A_GIMME Message Handling**: Variable-length list input for amplitude arrays
3. **State Management**: Coordinated control between different parameter methods
4. **Performance Optimization**: Wavetable synthesis with conditional processing
5. **Memory Management**: Dynamic allocation with proper cleanup
6. **Input Validation**: Robust handling of user input with helpful error messages

## See Also

- **Max Objects**: `cycle~`, `saw~`, `tri~` for basic waveforms
- **Related Externals**: `tide~` for LFO applications, `slewenv~` for envelope generation
- **Help File**: `harmosc~.maxhelp` for interactive examples
- **Development Notes**: `CLAUDE.md` for technical implementation details