---
layout: default
title: Setting up your computer
permalink: /setup/
---

# Setting up your computer

Most Isom Lab tools are desktop apps you start by double-clicking. Before the
first one will run, your computer needs one or two free pieces of software. **You
do this once**, not once per tool.

You do not need to know how to code. Apart from checking that an installer
worked, you will not need to type anything.

Mac and Windows are equally supported — every step below covers both. Each tool's
own guide says which of these steps it actually needs, and if it contradicts this
page, follow the tool: it knows about anything unusual it requires.

---

## 1. Miniforge

Most lab tools need **conda**, which gives each tool a private copy of Python and
everything it depends on, without disturbing anything else on your computer.
Miniforge is the free installer that provides it.

A few tools need none of this. Their guides say so, and send you to
[Tools that don't need conda](#tools-that-dont-need-conda) instead.

### Mac

1. Go to **[conda-forge.org/download](https://conda-forge.org/download/)**.
2. Download the macOS **`.pkg`** installer for your Mac: **Apple Silicon**
   (`arm64`) for any Mac from about 2020 onwards, **Intel** (`x86_64`) for older
   ones. Not sure which you have? Apple menu → **About This Mac**; a "Chip" line
   reading M1, M2, M3 or M4 means Apple Silicon.
3. Double-click the `.pkg` and click through, **accepting the defaults**.
4. Check it worked: quit Terminal if it's open, then open a fresh one (⌘-Space,
   type `Terminal`, Enter). You should see **`(base)`** at the start of the line.
   Close it again.

### Windows

1. Go to **[conda-forge.org/download](https://conda-forge.org/download/)**.
2. Under **Windows**, download **`Miniforge3-Windows-x86_64.exe`** — the right
   choice for every normal PC, including new ones.
3. Double-click it. If Windows says "Windows protected your PC", that is
   SmartScreen reacting to a freshly downloaded program, not a sign of trouble:
   click **More info ▸ Run anyway**.
4. Click through, **accepting the defaults**. Two screens are worth a glance:
   - **Install for: Just Me (recommended)** — leave it. It installs into your own
     user folder and needs no administrator password, which is what makes this
     possible on a shared or locked-down instrument PC.
   - The checkboxes near the end, about PATH and the default Python — leave them
     alone. The lab's launchers find conda by themselves.
5. Check it worked: open the Start menu, type `Miniforge`, and open **Miniforge
   Prompt**. A window appears with **`(base)`** at the start of the line. That is
   all it is for — close it again.

---

## 2. GitHub Desktop (optional)

Every tool lives in a repository on GitHub. You can always download one as a ZIP
straight from its page, and each tool's guide shows you how. GitHub Desktop does
the same job with buttons, and updates the tool later the same way.

1. Go to **[desktop.github.com](https://desktop.github.com)**, download, install.
2. Open it. You only need to sign in for the private tools described below.

On Windows there is one extra advantage: files GitHub Desktop downloads are never
"blocked", so you will not meet SmartScreen again.

> **You do not need Git.** Git is a separate command-line tool. None of the lab's
> apps ask you to use it, and GitHub Desktop already includes what it needs.

---

## Private tools

Most lab tools are public — anyone can download them, no account required. A few
are private, and those need a free GitHub account that has been given access.

1. Go to **[github.com](https://github.com)** and click **Sign up**. Enter your
   email, a password, and a **username** — keep that handy. Finish signing up and
   verify your email.
2. Send your GitHub username to Dan at
   <a href="&#109;&#97;&#105;&#108;&#116;&#111;&#58;&#100;&#105;&#115;&#111;&#109;&#64;&#109;&#105;&#97;&#109;&#105;&#46;&#101;&#100;&#117;">&#100;&#105;&#115;&#111;&#109;<span>&#64;</span>&#109;&#105;&#97;&#109;&#105;<span>&#46;</span>&#101;&#100;&#117;</a>.
   He adds you to the repository, and you get an email titled "[isomlab] invited
   you to collaborate". Open it and click **Accept invitation**.

You cannot download a private tool until that invitation is accepted.

---

## Tools that don't need conda

A few lab tools are written in nothing but Python's own standard library. They
have no environment to build; they just need a Python 3 on the computer.

- **Mac:** you already have one. Nothing to do.
- **Windows:** install Python from
  **[python.org/downloads](https://www.python.org/downloads/)** and tick **"Add
  Python to PATH"** on the first screen of the installer.

Both of those include **Tk**, the toolkit these tools use to open their windows,
so there is nothing further to install. If you installed Miniforge for another
tool, its Python includes Tk as well.

---

## You're set up

That is the shared groundwork. From here, follow the individual tool's own
getting-started guide, which is usually: download the folder, double-click the
launcher inside `launchers/`, and wait a few minutes the first time while it
builds its own private environment.

The tools are listed on the **[main page]({{ site.baseurl }}/)**.

---

## Trouble?

- **After installing Miniforge you don't see `(base)`** — close the window and open
  a brand-new one; the installer only affects windows opened afterwards. On Mac, if
  it still doesn't appear, run `source ~/miniforge3/bin/activate` once. On Windows,
  make sure you opened **Miniforge Prompt** from the Start menu, and not the
  ordinary Command Prompt or PowerShell — those don't know about conda.
- **The Miniforge installer asks for an administrator password (Windows)** — you
  chose "All Users" rather than **Just Me**. Go back and pick Just Me; it needs no
  password and works exactly the same.
- **"Windows protected your PC"** — SmartScreen, on anything downloaded from the
  internet. **More info ▸ Run anyway**. To avoid it on a downloaded ZIP,
  right-click the ZIP → **Properties** → tick **Unblock** → OK **before**
  extracting.
- **A launcher window flashes open and shuts (Windows)** — something failed too
  fast to read. Right-click the launcher → **Copy as path**, open **Miniforge
  Prompt**, paste, press Enter. The same thing runs, but the message stays on
  screen. Send what it says to Dan at
  <a href="&#109;&#97;&#105;&#108;&#116;&#111;&#58;&#100;&#105;&#115;&#111;&#109;&#64;&#109;&#105;&#97;&#109;&#105;&#46;&#101;&#100;&#117;">&#100;&#105;&#115;&#111;&#109;<span>&#64;</span>&#109;&#105;&#97;&#109;&#105;<span>&#46;</span>&#101;&#100;&#117;</a>.
- **You already have Anaconda or Miniconda** — that works too; you don't need
  Miniforge as well.
- **Didn't get the invitation email** — check spam, and confirm Dan used the exact
  username you signed up with. Pending invitations also appear at
  [github.com/notifications](https://github.com/notifications).
