---
layout: default
title: Setting up your computer
permalink: /setup/
---

# Setting up your computer

Most Isom Lab tools are desktop apps you start by double-clicking. Before the
first one will run, your computer needs a few free pieces of software. **You do
this once**, not once per tool.

You do **not** need to know how to code, and — apart from checking that an
installer worked — you will not need to type anything.

> **Mac and Windows are equally supported.** Every step below covers both. If a
> tool's own guide contradicts this page, follow the tool's guide: it knows about
> anything unusual it needs.

Each tool's own `docs/INSTALL.md` tells you which of these steps it actually
needs, and adds anything specific to it.

---

## 1. A GitHub account *(only for private tools)*

Some lab tools are in private repositories. To download one, you need a free
GitHub account that has been given access.

1. Go to **[github.com](https://github.com)** and click **Sign up**. Enter your
   email, a password, and a **username** — keep that handy, you'll send it to Dan.
   Finish signing up and **verify your email**.
2. **Send your GitHub username to Dan.** He adds you to the repository, and you'll
   get an email titled *"[isomlab] invited you to collaborate"*. Open it and click
   **Accept invitation**.

You can't download a private tool until that invitation is accepted. Tools that
are public need none of this.

---

## 2. Miniforge *(for tools that say they need conda)*

Miniforge gives you a private copy of Python and everything a tool needs, without
disturbing anything else on your computer. Some lab tools are pure Python and
need nothing at all — their guides say so. Skip this step for those.

### Mac

1. Go to **[conda-forge.org/download](https://conda-forge.org/download/)**.
2. Choose the **macOS** installer that matches your Mac:
   - **Apple Silicon** (`arm64`) — Macs from about 2020 onwards. *(Not sure? Apple
     menu → **About This Mac**. "Chip: Apple M1/M2/M3/M4" means Apple Silicon.)*
   - **Intel** (`x86_64`) — older Macs, listed as "Intel Core i5/i7".
   - Take the **`.pkg`** installer: that's the click-through one.
3. Double-click the downloaded `.pkg` and click through, **accepting the defaults**.
4. Check it worked: **quit Terminal if it's open, then open a fresh one**
   (⌘-Space, type `Terminal`, Enter). You should see **`(base)`** at the start of
   the line. Close it again.

### Windows

1. Go to **[conda-forge.org/download](https://conda-forge.org/download/)**.
2. Under **Windows**, download **`Miniforge3-Windows-x86_64.exe`** — the right
   choice for every normal PC, including new ones.
3. Double-click it. Windows may say *"Windows protected your PC"*; that's
   SmartScreen reacting to a freshly downloaded program, not a sign of trouble.
   Click **More info ▸ Run anyway**.
4. Click through, **accepting the defaults**. Two screens are worth a glance:
   - **"Install for: Just Me (recommended)"** — leave it. It installs into your own
     user folder and **needs no administrator password**, which is what makes this
     possible on a shared or locked-down instrument PC.
   - The checkboxes near the end, about PATH and the default Python — **leave them
     alone**. The lab's launchers find conda by themselves.
5. Check it worked: open the Start menu, type `Miniforge`, and open **Miniforge
   Prompt**. A window appears with **`(base)`** at the start of the line. That is
   all it's for — close it again. **You will not need to type anything.**

---

## 3. GitHub Desktop *(optional, both platforms)*

This downloads a tool, and later updates it, with buttons instead of typing.

1. Go to **[desktop.github.com](https://desktop.github.com)**, download, install.
2. Open it and **sign in** with your GitHub account.

You can skip this and download a ZIP from the repository page instead — each
tool's guide shows both. One advantage of GitHub Desktop on Windows: files it
downloads are never "blocked" by Windows, so you won't meet SmartScreen again.

> **You do not need Git.** Git is a separate command-line tool. None of the lab's
> apps ask you to use it, and GitHub Desktop already includes what it needs.

---

## You're set up

That's the shared groundwork. From here, follow the individual tool's
`docs/getting_started.md`, which is usually: download the folder, double-click the
launcher inside `launchers/`, and wait a few minutes the first time while it
builds its own private environment.

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
- **Which Mac installer?** Apple menu → About This Mac → the "Chip"/"Processor"
  line tells you Apple Silicon or Intel.
- **"Windows protected your PC"** — SmartScreen, on anything downloaded from the
  internet. **More info ▸ Run anyway**. To avoid it on a downloaded ZIP,
  right-click the ZIP → **Properties** → tick **Unblock** → OK *before* extracting.
- **A launcher window flashes open and shuts (Windows)** — something failed too
  fast to read. Right-click the launcher → **Copy as path**, open **Miniforge
  Prompt**, paste, press Enter. The same thing runs, but the message stays on
  screen. Send Dan what it says.
- **You already have Anaconda or Miniconda** — that works too; you don't need
  Miniforge as well.
- **Didn't get the invitation email** — check spam, and confirm Dan used the exact
  username you signed up with. Pending invitations also appear at
  [github.com/notifications](https://github.com/notifications).
