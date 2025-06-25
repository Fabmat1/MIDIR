import json
import os
import tkinter as tk
import tkinter.filedialog as fd
import tkinter.ttk as ttk
from PIL import ImageTk, Image
import sv_ttk

from src.GUI import DataReductionGUI

if __name__ == "__main__":
	root = DataReductionGUI()
	sv_ttk.set_theme("light")
	root.mainloop()