# MIDIR – Multi-Instrument Data Input Reduce-inator

**MIDIR** is a flexible and modular data reduction toolkit designed to streamline the processing of spectral data from long slit spectrographs. With a focus on usability, reproducibility, and performance, MIDIR provides a graphical user interface and robust backend for handling a variety of input formats and reduction tasks. the goal of MIDIR is to provide a framework that can effortlessly be expanded to work with any long slit spectrograph as a reduction tool.

## Features

- Graphical interface for configuring reduction parameters
- Modular design for adding support for new instruments
- Support for cropping, bias correction, flat-fielding, and more
- Customizable presets for recurring workflows
- Real-time progress tracking

## Installation

### 1. Clone the repository:
```bash
git clone https://github.com/Fabmat1/MIDIR.git
cd MIDIR
```
### 2. Install the required Python dependencies:

```bash
pip install -r requirements.txt
```

### 3. Compile required components:

```bash
make
```

> ⚠️ You must run make before launching MIDIR to ensure all components are correctly built.

## Usage

To start the application:

```bash
python main.py
```

The GUI will guide you through selecting your input files, configuring reduction options, and saving or loading presets.

[Read the full instructions on using the GUI here](INSTRUCTIONS.md)

## Folder Structure

```bash
MIDIR/
├── main.py                  # Entry point for the application
├── requirements.txt         # Python dependencies
├── Makefile                 # Build script for native components
├── src/                     # Source code
│   ├── options.py           # Reduction configuration
│   ├── frames.py            # Frame processing logic
│   └── ...                  # Additional modules
└── README.md                # Project documentation
```

## Contributing

Contributions are welcome! Please open an issue or pull request if you would like to suggest improvements, report bugs, or add features.

## License

This project is licensed under the MIT License. See LICENSE for details.
