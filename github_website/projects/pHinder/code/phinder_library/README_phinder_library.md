# phinder_library

This folder contains the core Python modules used by the Isom Lab’s [pHinder](../../) program. These scripts implement the main computational functionality and are required by all user-facing components, including `phinder_gui` and `phinder_command_line`.

---

## 🚀 How to Use

To use this code in your own environment, follow the steps below:

### Step 1: Download the code

Download this folder manually, or clone the full repository:

```
git clone https://github.com/dangerisom/IsomLabPrivate.git
```

Then navigate to the subfolder containing the core library:

```
cd IsomLabPrivate/projects/pHinder/code/phinder_library
```

---

### Step 2: Add the folder to your Python path

To use this code as a dependency in your own scripts or in the GUI, you must add this folder to your Python path.

**Option A – Add temporarily inside your script:**

```python
import sys
sys.path.append('/path/to/phinder_library')
```

**Option B – Add permanently via environment variable:**

**On macOS/Linux:**

```
export PYTHONPATH=$PYTHONPATH:/path/to/phinder_library
```

**On Windows:**

```
set PYTHONPATH=%PYTHONPATH%;C:\path\to\phinder_library
```

Replace `/path/to/phinder_library` with the full path on your system.

---

### Step 3: Install required third-party modules listed in `3rd_party_software.txt`

All external dependencies used by this codebase are listed in `3rd_party_software.txt`.

Install the packages using either pip or conda:

**Using pip:**

```
pip install numpy matplotlib
```

**Using conda:**

```
conda install numpy matplotlib
```

> 📎 Check the contents of `3rd_party_software.txt` for the exact list of packages.

---

## 📁 Folder Contents

- Core Python `.py` source files used by the pHinder suite  
- `3rd_party_software.txt` – List of required third-party packages  
- `README.md` – This setup and usage guide  

---

## 📬 Questions?

If you have questions about using this library or its dependencies, please contact the Isom Lab or open an issue in the repository.
