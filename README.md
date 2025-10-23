# McGill Rocket Team Steady-Unsteady-Flight-Model

This repository contains the computational models (Steady-State and Unsteady) for the nitrous oxide/paraffin wax hybrid rocket engine.

## 🚀 Quick Start Guide

Follow these steps to get the code running on your machine.

### 1. Prerequisites

*   **Python 3.13**
*   **Git**

**Check your Python version:**
```bash
python --version
```

### 2. Clone this repository
```bash
git clone https://github.com/Leo11235/Steady-Unsteady-Flight-Model
cd Steady-Unsteady-Flight-Model
```

### 3. Set up your environment
Windows
```bash
# Create virtual environment
python -m venv venv

# Activate it
venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

Mac/Linux
```bash
# Create virtual environment
python3 -m venv venv

# Activate it  
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### 4. Verify your setup
Run the verification test to make sure everything works:
```bash
python backend/steady_venv_testing.py
```
You should see simulation loops running with apogee calculations, ending with:
```bash
LOOP 9 - FINAL APOGEE: 13715.2 meters (0.81256 meters off from target)
```
If you see this, your setup is complete and working!
