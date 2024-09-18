# ComplexWaveletFilter
The complex wavelet filtering approach applied phasor transformed TCSPC FLIM data used in Fahim and Marcus et al. JCB 2024.

## Installation
To set up the environment and run the scripts, follow these steps:

### 0.0. Installation Requirements
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
Once the dependencies are installed, you can run the main script:

```bash
python ComplexWaveletFilter.py
```

A window will popup to select the directory containing .tif files with calibrated and unfiltered phasor coordiantes (G.tif and S.tif) as well as an intensity image (intensity.tif). A sample dataset is provided. If using the sample dataset, unzip sample_data.zip and select the directory sample_data.

Next, enter the harmonic used for the phasor transformation

Next, enter the expected lifetime of the fluoraphore being imaged

Last, enter the desired levels of filtering

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
