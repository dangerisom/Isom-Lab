# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
from datetime import datetime
import cv2
import numpy as np
import csv
import sys
from typing import List, Tuple, Optional

# Default is usually 1000
sys.setrecursionlimit(10000)
print(sys.getrecursionlimit())

class Data:
    def __init__(self):
        self.label = ""
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.mixed = False

class ImagePixelOverlapApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Vesicle Co-localization Overlap GUI")
        # Resize window so all fields fit comfortably
        self.geometry("900x700")
        self.minsize(900, 700)

        # List of paths to coordinate csv files
        self.csv_paths = []

        # List of paths to calculated contour image files
        self.contour_image_paths = set()

        # Path to contour image overlay
        self.contour_image_overlay_path = ""

        # Base directories
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.default_dir = self.script_dir

        # Threshold and contour parameters
        self.t1_min = 50
        self.t1_int = 150
        self.t1_max = 255
        self.t2_min = 50
        self.t2_int = 150
        self.t2_max = 255
        self.contour_width = 1
        self.contour_area_min = 1
        self.contour_area_max = 1000
        self.contour_color_channel_1 = (0, 0, 255) # red
        self.contour_color_channel_2 = (0, 0, 255) # red

        # Co-localization and overlap results
        self.total1 = 0
        self.total2 = 0
        self.overlap_count = 0
        self.frac1 = 0.0
        self.frac2 = 0.0

        # === Image Channel 1 ===
        self.frame_channel1 = ttk.LabelFrame(self, text="Image Channel 1")
        self.frame_channel1.pack(padx=20, pady=(20, 10), fill="x")
        self.btn_channel1 = ttk.Button(
            self.frame_channel1,
            text="Load Channel 1",
            command=self.load_channel1
        )
        self.btn_channel1.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        self.lbl_channel1 = ttk.Label(
            self.frame_channel1,
            text="No file selected",
            foreground="gray"
        )
        self.lbl_channel1.grid(row=0, column=1, padx=10, pady=10, sticky="w")

        # === Image Channel 2 ===
        self.frame_channel2 = ttk.LabelFrame(self, text="Image Channel 2")
        self.frame_channel2.pack(padx=20, pady=(0, 10), fill="x")
        self.btn_channel2 = ttk.Button(
            self.frame_channel2,
            text="Load Channel 2",
            command=self.load_channel2
        )
        self.btn_channel2.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        self.lbl_channel2 = ttk.Label(
            self.frame_channel2,
            text="No file selected",
            foreground="gray"
        )
        self.lbl_channel2.grid(row=0, column=1, padx=10, pady=10, sticky="w")

        # === Save Output ===
        self.frame_save = ttk.LabelFrame(self, text="Save Output Directory")
        self.frame_save.pack(padx=20, pady=(0, 20), fill="x")
        self.btn_save = ttk.Button(
            self.frame_save,
            text="Select Save Location",
            command=self.select_save_directory
        )
        self.btn_save.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        self.lbl_save = ttk.Label(
            self.frame_save,
            text="No directory selected",
            foreground="gray"
        )
        self.lbl_save.grid(row=0, column=1, padx=10, pady=10, sticky="w")

        # ─── Thresholds Frame ───────────────────────────────────────────────
        self.frame_thresholds = ttk.LabelFrame(self, text="Channel Thresholds")
        self.frame_thresholds.pack(padx=20, pady=(0, 20), fill="x")

        # Channel 1 Min / Int / Max
        ttk.Label(self.frame_thresholds, text="Ch1 Min:")\
            .grid(row=0, column=0, sticky="w", padx=10, pady=5)
        self.entry_t1_min = ttk.Entry(self.frame_thresholds, width=6)
        self.entry_t1_min.insert(0, str(self.t1_min))
        self.entry_t1_min.grid(row=0, column=1, sticky="w")

        ttk.Label(self.frame_thresholds, text="Ch1 Int:")\
            .grid(row=0, column=2, sticky="w", padx=10)
        self.entry_t1_int = ttk.Entry(self.frame_thresholds, width=6)
        self.entry_t1_int.insert(0, str(self.t1_int))
        self.entry_t1_int.grid(row=0, column=3, sticky="w")

        ttk.Label(self.frame_thresholds, text="Ch1 Max:")\
            .grid(row=0, column=4, sticky="w", padx=10)
        self.entry_t1_max = ttk.Entry(self.frame_thresholds, width=6)
        self.entry_t1_max.insert(0, str(self.t1_max))
        self.entry_t1_max.grid(row=0, column=5, sticky="w")

        # Channel 2 Min / Int / Max
        ttk.Label(self.frame_thresholds, text="Ch2 Min:")\
            .grid(row=1, column=0, sticky="w", padx=10, pady=5)
        self.entry_t2_min = ttk.Entry(self.frame_thresholds, width=6)
        self.entry_t2_min.insert(0, str(self.t2_min))
        self.entry_t2_min.grid(row=1, column=1, sticky="w")

        ttk.Label(self.frame_thresholds, text="Ch2 Int:")\
            .grid(row=1, column=2, sticky="w", padx=10)
        self.entry_t2_int = ttk.Entry(self.frame_thresholds, width=6)
        self.entry_t2_int.insert(0, str(self.t2_int))
        self.entry_t2_int.grid(row=1, column=3, sticky="w")

        ttk.Label(self.frame_thresholds, text="Ch2 Max:")\
            .grid(row=1, column=4, sticky="w", padx=10)
        self.entry_t2_max = ttk.Entry(self.frame_thresholds, width=6)
        self.entry_t2_max.insert(0, str(self.t2_max))
        self.entry_t2_max.grid(row=1, column=5, sticky="w")

        # Preview button spanning both rows
        self.btn_preview = ttk.Button(
            self.frame_thresholds,
            text="Preview Threshold",
            command=self.on_preview_threshold
        )
        self.btn_preview.grid(row=0, column=6, rowspan=2, padx=10, pady=5, sticky="nsew")

        
        # === Threshold Settings ===
        self.frame_params = ttk.LabelFrame(self, text="Contour Settings")
        self.frame_params.pack(padx=20, pady=(0, 10), fill="x")

        # Row 1: Contour Width
        ttk.Label(self.frame_params, text="Contour Width:")\
            .grid(row=1, column=0, padx=10, pady=5, sticky="w")
        self.entry_contour_width = ttk.Entry(self.frame_params, width=8)
        self.entry_contour_width.insert(0, str(self.contour_width))
        self.entry_contour_width.grid(row=1, column=1, padx=10, pady=5, sticky="w")

        # Row 2: Area Min / Max
        ttk.Label(self.frame_params, text="Area Min:")\
            .grid(row=2, column=0, padx=10, pady=5, sticky="w")
        self.entry_area_min = ttk.Entry(self.frame_params, width=8)
        self.entry_area_min.insert(0, str(self.contour_area_min))
        self.entry_area_min.grid(row=2, column=1, padx=10, pady=5, sticky="w")

        ttk.Label(self.frame_params, text="Area Max:")\
            .grid(row=2, column=2, padx=10, pady=5, sticky="w")
        self.entry_area_max = ttk.Entry(self.frame_params, width=8)
        self.entry_area_max.insert(0, str(self.contour_area_max))
        self.entry_area_max.grid(row=2, column=3, padx=10, pady=5, sticky="w")

        # Row 4: Process Button (left‑justified, no spanning)
        self.btn_threshold_and_contours = ttk.Button(
            self.frame_params,
            text="Calculate Pixel Overlap",
            command=self.on_threshold_and_contours
            # command=self.on_threshold
        )
        self.btn_threshold_and_contours.grid(
            row=4, column=0, sticky="w", padx=10, pady=(10,5)
        )

        # Internal variables
        self.channel1_path = None
        self.channel2_path = None
        self.save_dir = None
        self.vertices = []

    def load_channel1(self):
        path = filedialog.askopenfilename(
            title="Select Image Channel 1",
            filetypes=[("Image files", "*.png *.jpg *.tif *.tiff *.bmp")],
            initialdir=self.default_dir
        )
        if path:
            self.channel1_path = path
            self.lbl_channel1.config(text=os.path.basename(path), foreground="green")

    def load_channel2(self):
        path = filedialog.askopenfilename(
            title="Select Image Channel 2",
            filetypes=[("Image files", "*.png *.jpg *.tif *.tiff *.bmp")],
            initialdir=self.default_dir
        )
        if path:
            self.channel2_path = path
            self.lbl_channel2.config(text=os.path.basename(path), foreground="green")

    def select_save_directory(self):
        directory = filedialog.askdirectory(initialdir=self.default_dir)
        if directory:
            timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            save_path = os.path.join(directory, f"output_{timestamp}")
            os.makedirs(save_path, exist_ok=True)
            self.save_dir = save_path
            self.lbl_save.config(text=self.save_dir, foreground="green")


    def on_preview_threshold(self):
        """
        Preview thresholding for both channels using current GUI values.
        Shows OpenCV windows; does not save files or alter state.
        """
        if not self.channel1_path or not self.channel2_path:
            messagebox.showwarning("Missing Channels", "Please load both image channels before previewing.")
            return

        # use a directory even if not saving (required by signature)
        out_dir = self.save_dir if getattr(self, "save_dir", None) else os.getcwd()

        try:
            # Read contour params
            self.contour_width        = int(self.entry_contour_width.get())
            self.contour_area_min     = int(self.entry_area_min.get())
            self.contour_area_max     = int(self.entry_area_max.get())

            # Colors
            self.contour_color_channel_1 = (0, 0, 225) #parse_bgr(self.entry_contour_color_channel_1.get())
            self.contour_color_channel_2 = (0, 0, 225) #parse_bgr(self.entry_contour_color_channel_2.get())

            # Thresholds
            t1_min = int(self.entry_t1_min.get()); t1_int = int(self.entry_t1_int.get()); t1_max = int(self.entry_t1_max.get())
            t2_min = int(self.entry_t2_min.get()); t2_int = int(self.entry_t2_int.get()); t2_max = int(self.entry_t2_max.get())
        except Exception as e:
            messagebox.showerror("Invalid Input", f"Please enter valid numeric values and BGR tuples.\n{e}")
            return

        # Run previews (display=True, save=False)
        self.threshold_and_contours(self.channel1_path, self.contour_color_channel_1, out_dir,
                                    t1_min, t1_int, t1_max, return_mask=False, display=True, save=False)
        self.threshold_and_contours(self.channel2_path, self.contour_color_channel_2, out_dir,
                                    t2_min, t2_int, t2_max, return_mask=False, display=True, save=False)
        
    def on_threshold_and_contours(self):

        if not self.channel1_path or not self.channel2_path:
            messagebox.showwarning(
                "Missing Channels",
                "Please load both image channels before processing."
            )
            return
        if not self.save_dir:
            messagebox.showwarning(
                "No Save Directory",
                "Please select a save directory before processing."
            )
            return

        # ─── Reset any previous outputs ──────────────────────────────────────
        self.csv_paths.clear()
        self.contour_image_paths.clear()
        self.contour_image_overlay_path = ""

        # ─── Read values from the GUI ────────────────────────────────────────

        # Read contour params as before…
        self.contour_width    = int(self.entry_contour_width.get())
        self.contour_area_min = int(self.entry_area_min.get())
        self.contour_area_max = int(self.entry_area_max.get())

        # Read _independent_ thresholds:
        try:
            t1_min = self.t1_min = int(self.entry_t1_min.get())
            t1_int = self.t1_int = int(self.entry_t1_int.get())
            t1_max = self.t1_max = int(self.entry_t1_max.get())

            t2_min = self.t2_min = int(self.entry_t2_min.get())
            t2_int = self.t2_int = int(self.entry_t2_int.get())
            t2_max = self.t2_max =int(self.entry_t2_max.get())
        except ValueError as e:
            messagebox.showerror("Invalid Threshold",
                                 f"Enter valid integers:\n{e}")
            return

        # Read contour parameters
        try:
            self.contour_width        = int(self.entry_contour_width.get())
            self.contour_area_min     = int(self.entry_area_min.get())
            self.contour_area_max     = int(self.entry_area_max.get())

            # def parse_bgr(s):
            #     parts = s.strip("()[] ").split(",")
            #     if len(parts) != 3:
            #         raise ValueError(f"Expected 3 values but got {len(parts)}")
            #     return tuple(int(p) for p in parts)

            self.contour_color_channel_1 = (0, 0, 255) #parse_bgr(self.entry_contour_color_channel_1.get())
            self.contour_color_channel_2 = (0, 0, 255) #parse_bgr(self.entry_contour_color_channel_2.get())

        except ValueError as e:
            messagebox.showerror(
                "Invalid Input",
                f"Please enter valid numeric values or BGR color tuples:\n{e}"
            )
            return

        # Generate masks for each channel
        mask1 = self.threshold_and_contours(
            self.channel1_path,
            self.contour_color_channel_1,
            self.save_dir,
            t1_min, t1_int, t1_max,
            return_mask=True
        )
        mask2 = self.threshold_and_contours(
            self.channel2_path,
            self.contour_color_channel_2,
            self.save_dir,
            t2_min, t2_int, t2_max,
            return_mask=True
        )

        # 1. Generate binary masks for each channel
        mask1 = self.threshold_and_contours(
            self.channel1_path,
            self.contour_color_channel_1,
            self.save_dir,
            t1_min, t1_int, t1_max,
            return_mask=True
        )
        mask2 = self.threshold_and_contours(
            self.channel2_path,
            self.contour_color_channel_2,
            self.save_dir,
            t2_min, t2_int, t2_max,
            return_mask=True
        )

        # 2. Compute overlap
        overlap_mask = cv2.bitwise_and(mask1, mask2)

        # 3–5. Counts & fractions
        self.total1       = int(np.count_nonzero(mask1))
        self.total2       = int(np.count_nonzero(mask2))
        self.overlap_count= int(np.count_nonzero(overlap_mask))
        self.frac1        = (self.overlap_count/self.total1) if self.total1 else 0
        self.frac2        = (self.overlap_count/self.total2) if self.total2 else 0

        # 6. Save “overlap on masked” for each channel
        for path, mask, color in [
            (self.channel1_path, mask1, self.contour_color_channel_1),
            (self.channel2_path, mask2, self.contour_color_channel_2)
        ]:
            base = os.path.splitext(os.path.basename(path))[0]
            # reload original color image
            img = cv2.imread(path)
            # recreate masked image
            masked = cv2.bitwise_and(img, img, mask=mask)
            # overlay overlap pixels
            overlay = masked.copy()
            overlay[overlap_mask > 0] = color
            # save
            fn = os.path.join(self.save_dir, f"{base}_overlap_on_masked.png")
            cv2.imwrite(fn, overlay)

        # 7. Finally, show results
        self.record_parameters_and_results()
        messagebox.showinfo(
            "Overlap Results",
            (f"Overlap pixels: {self.overlap_count}\n"
             f"Total Ch1 pixels: {self.total1}\n"
             f"Total Ch2 pixels: {self.total2}\n"
             f"Frac Ch1: {self.frac1:.3f}\n"
             f"Frac Ch2: {self.frac2:.3f}")
        )

    def threshold_and_contours(self, path, contour_color, output_dir,
                               t_min, t_int, t_max,
                               return_mask=False, display=False, save=True):
        """
        If return_mask=True, returns the binary mask (uint8) of all selected pixels.
        Otherwise behaves exactly as before.
        """
        img_gray  = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
        img_color = cv2.imread(path)
        if img_gray is None or img_color is None:
            raise ValueError(f"Failed to read image: {path}")

        base      = os.path.splitext(os.path.basename(path))[0]
        overlay   = img_color.copy()
        cont_only = np.zeros_like(img_color)

        # build two submasks, then combine
        mask_low  = np.zeros_like(img_gray, dtype=np.uint8)
        mask_high = np.zeros_like(img_gray, dtype=np.uint8)

        for lo, hi, mask in [(t_min, t_int, mask_low),
                             (t_int, t_max, mask_high)]:

            _, thresh = cv2.threshold(img_gray, lo, hi, cv2.THRESH_BINARY)
            contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL,
                                           cv2.CHAIN_APPROX_SIMPLE)

            for cnt in contours:
                area = cv2.contourArea(cnt)
                if not (self.contour_area_min < area < self.contour_area_max):
                    continue
                # outlines
                cv2.drawContours(overlay,    [cnt], -1,
                                 contour_color, self.contour_width)
                cv2.drawContours(cont_only,  [cnt], -1,
                                 contour_color, self.contour_width)
                # fill mask
                cv2.drawContours(mask, [cnt], -1, 255, thickness=-1)

        # final binary mask of everything selected
        mask_all = cv2.bitwise_or(mask_low, mask_high)

        # save masked PNG
        masked = cv2.bitwise_and(img_color, img_color, mask=mask_all)
        if save:
            cv2.imwrite(os.path.join(output_dir, f"{base}_masked.png"), masked)

        # save overlay + contour-only as before
        if save:
            cv2.imwrite(os.path.join(output_dir, f"{base}_overlay.png"), overlay)
        if save:
            cv2.imwrite(os.path.join(output_dir, f"{base}_contours.png"), cont_only)
        if save:
            self.contour_image_paths.add(f"{base}_contours.png")

        if display:
            # cv2.imshow(f"{base} Masked", masked)
            cv2.imshow(f"{base} Overlay", overlay)
            # cv2.imshow(f"{base} Contours", cont_only)
            cv2.waitKey(0)
            cv2.destroyAllWindows()

        if return_mask:
            return mask_all

    def record_parameters_and_results(self):
        from datetime import datetime as _dt
        import os

        # Timestamp & log file
        now = _dt.now()
        timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
        log_name = f"vesicle_coloc_overlap_log_{now.strftime('%Y%m%d_%H%M%S')}.txt"
        log_path = os.path.join(self.save_dir, log_name)

        # Derive base filenames for channel outputs
        base1 = os.path.splitext(os.path.basename(self.channel1_path))[0]
        base2 = os.path.splitext(os.path.basename(self.channel2_path))[0]

        files = {
            "Ch1 Masked":           f"{base1}_masked.png",
            "Ch1 Overlap-Masked":   f"{base1}_overlap_on_masked.png",
            "Ch1 Overlay":          f"{base1}_overlay.png",
            "Ch1 Contours":         f"{base1}_contours.png",
            "Ch2 Masked":           f"{base2}_masked.png",
            "Ch2 Overlap-Masked":   f"{base2}_overlap_on_masked.png",
            "Ch2 Overlay":          f"{base2}_overlay.png",
            "Ch2 Contours":         f"{base2}_contours.png",
            "Log File":             log_name
        }

        # Collect parameters
        params = {
            "Channel 1 Min":           getattr(self, "t1_min", None),
            "Channel 1 Intermediate":  getattr(self, "t1_int", None),
            "Channel 1 Max":           getattr(self, "t1_max", None),
            "Channel 2 Min":           getattr(self, "t2_min", None),
            "Channel 2 Intermediate":  getattr(self, "t2_int", None),
            "Channel 2 Max":           getattr(self, "t2_max", None),
            "Contour Width":           self.contour_width,
            "Contour Area Min":        self.contour_area_min,
            "Contour Area Max":        self.contour_area_max
        }

        # Results summary
        results = {
            "Total Pixels Image 1":    self.total1,
            "Total Pixels Image 2":    self.total2,
            "Overlapping Pixels":      self.overlap_count,
            "Fraction Overlap Ch1":    self.frac1,
            "Fraction Overlap Ch2":    self.frac2
        }

        # Build log lines
        lines = [
            "=== Coloc Pixel Overlay Log ===",
            f"Date: {timestamp}", "",
            "Inputs:",
            f"  Image 1: {self.channel1_path}",
            f"  Image 2: {self.channel2_path}", ""
        ]

        lines.append("Parameters:")
        for k, v in params.items():
            lines.append(f"  {k}: {v}")
        lines.append("")

        lines.append("Results:")
        lines.append(f"  Total Pixels Image 1: {results['Total Pixels Image 1']}")
        lines.append(f"  Total Pixels Image 2: {results['Total Pixels Image 2']}")
        lines.append(f"  Overlapping Pixels:   {results['Overlapping Pixels']}")
        lines.append(f"  Fraction Overlap Ch1: {results['Fraction Overlap Ch1']:.4f}")
        lines.append(f"  Fraction Overlap Ch2: {results['Fraction Overlap Ch2']:.4f}")
        lines.append("")

        lines.append("Output Files:")
        for desc, fname in files.items():
            # include full path for everything except log_name (we include log_path below)
            path = log_path if desc == "Log File" else os.path.join(self.save_dir, fname)
            lines.append(f"  {desc}: {path}")

        # Write out
        with open(log_path, "w") as lf:
            lf.write("\n".join(lines))

        # Also print to console
        print("\n".join(lines))
        print(f"[INFO] Wrote parameter/results log to {log_path}")

if __name__ == "__main__":
    app = ImagePixelOverlapApp()
    app.mainloop()
