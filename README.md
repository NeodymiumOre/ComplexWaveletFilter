# ComplexWaveletFilter
The complex wavelet filtering approach applied to phasor transformed TCSPC FLIM data used in Fahim and Marcus et al. JCB 2025. (https://doi.org/10.1083/jcb.202311105)

## Installation
To set up the environment and run the scripts, follow these steps:

### 0. Installation Requirements
Before you do anything, you’ll want to make sure your computer is properly set up for installation. This first installation version only works for Mac. Windows and Linux installation instructions will be included in later versions...

In order to properly install FLIMagePy, you will need the following:

- macOS 10.9 (Mavericks) or later
- Python 3.6 or later
- git 2.0 or later

Check the Python version installed on your computer by entering the following in your terminal app
```bash
python3 --version
```

Check the Python version installed on your computer by entering the following in your terminal app
```bash
git --version
```

For mac users, I recomend using Homebrew to install python and git. You can find the documentation and installation instructions at https://brew.sh/. Once Homebrew is installed on your machine, install the latest versions of python and git:
```bash
brew install python
brew install git
```

### 1. Clone the Repository
Open your terminal app and enter the following
```bash
cd ~
git clone https://github.com/marcusjoshm/ComplexWaveletFilter.git
cd ComplexWaveletFilter
```

### 2. Set Up a Virtual Environment
Create and activate a virtual environment:

```bash
# Navigate to your desired directory for the virtual environment
cd ~/ComplexWaveletFilter

# Create the virtual environment
python3 -m venv venv

# Activate the virtual environment
# On macOS/Linux:
source venv/bin/activate
```

### 3. Install the Required Dependencies
With the virtual environment activated, navigate to your project directory and install the dependencies:

```bash
cd ~/ComplexWaveletFilter
pip install -r requirements.txt
```

### 4. Run the Script
First, run the complex wavelet filtering script:

```bash
python3 ComplexWaveletFilter.py
```

A window will popup to select the directory containing .tif files with calibrated and unfiltered phasor coordiantes (G.tif and S.tif) as well as an intensity image (intensity.tif). A sample dataset is provided.If using the sample dataset, unzip sample_data.zip and select the directory sample_data.

Next, enter the harmonic used for the phasor transformation. The sample dataset provided was transformed using a harmonic of 1.

Next, enter the expected lifetime of the fluoraphore being imaged. mNeonGreen is the fluorophore used for the sample dataset, which has an expected lifetime of 3.1ns

Last, enter the desired levels of filtering. Wavelet filtered phasor plots from Fahim and Marcus et al. JCB 2024 use 9 filtering levels.

Next, run the condensed phase segmentation script:

```bash
python3 CondensedPhaseGMM.py
```

A window will popup to select a .npz file containing filtered phasor coordinates. If this was generated using ComplexWaveletFilter.py, the .npz generated from this program will be located in a subdirectory called datasets.

Next, enter the harmonic and expected lifetime of the fluoraphore (this should be the same values you entered for the ComplexWaveletFilter.py program).

Next will be prompts for segmentation perameters. You will be promted to enter a multiplication factor that defines the size of the segmentation ROI, a shift factor to move the center of the ROI if needed, and a radius for a filtering ROI (see methods section for detailed description). For the sample dataset provided, use a multiplication factor of 4, a shift factor of 0.625, and a radius of 0.1 for complete segmentation of the condensed phase cluster.

### Dependencies
The following Python packages are required and are listed in the `requirements.txt` file:

- numpy>=1.21.0
- Pillow>=8.0.0
- dtcwt>=0.12.0
- matplotlib>=3.4.0
- tifffile>=2021.7.2
- scikit-learn>=0.24.0
- scipy>=1.6.0

### License
This project is licensed under the MIT License - see the `LICENSE` file for details.
