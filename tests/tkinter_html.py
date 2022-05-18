from tkinter import *
from tkhtmlview import *

import codecs

root = Tk()
root.geometry('400x400')
root.title('HTML in tkinter')

with open("html_outputs/plot_dk.html", "r", encoding='utf-8') as f:
    html_text= f.read()
    
# create HTMLLabel
label = HTMLLabel(root)
label.set_html(html_text)
label.pack(pady=15, padx=15, fill=BOTH)

root.mainloop()