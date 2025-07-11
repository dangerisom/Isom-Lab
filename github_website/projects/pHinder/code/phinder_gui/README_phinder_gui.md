# phinder_gui

This folder contains Python code for the `phinder_gui` module, part of the Isom Lab’s [pHinder](../../) project. These scripts support a graphical interface for structural analysis workflows based in computational geometry.

---

## 🚀 How to Use

To use this code in your own environment, follow the steps below:

### Step 1: Download the code

Download this folder manually, or clone the full repository:

```
git clone https://github.com/dangerisom/IsomLabPrivate.git
```

Then navigate to the subfolder containing the code:

```
cd IsomLabPrivate/projects/pHinder/code/phinder_gui
```

---

### Step 2: Add the folder to your Python path

In addition to this folder, you must also download the `phinder_library` folder (located in the same repository) and ensure it is in your Python path. The scripts in `phinder_gui` depend on modules defined in `phinder_library`.


To allow Python to recognize and import the modules in this folder, you must add it to your Python path.

**Option A – Add temporarily inside your script:**

```python
import sys
sys.path.append('/path/to/phinder_gui')
```

**Option B – Add permanently via environment variable:**

**On macOS/Linux:**

```
export PYTHONPATH=$PYTHONPATH:/path/to/phinder_gui
```

**On Windows:**

```
set PYTHONPATH=%PYTHONPATH%;C:\path\to\phinder_gui
```

Replace `/path/to/phinder_gui` with the full path on your system.

---

### Step 3: Install required third-party modules listed in `3rd_party_software.txt`

All external dependencies needed by this module are listed in `3rd_party_software.txt`.

You can install these packages using either pip or conda:

**Using pip:**

```
pip install numpy matplotlib PyQt5
```

**Using conda:**

```
conda install numpy matplotlib pyqt
```

> 📎 Be sure to check the actual contents of `3rd_party_software.py` for the full list of required modules.

---

## 📁 Folder Contents

- Python `.py` source files for the GUI and supporting functions  
- `3rd_party_software.txt` – Lists required third-party packages  
- `README_phinder_gui.md` – This setup and usage guide  

---

## 📬 Questions?

If you encounter any issues using this module or installing dependencies, please contact the Isom Lab or open an issue in the repository.
