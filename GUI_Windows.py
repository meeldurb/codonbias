"""Author: Melanie van den Bosch
Jorn van der Ent
Script for having multiple windows frames in a GUI"""

from __future__ import division
import math 
import Tkinter as tk

class OligoDatabase(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        container = tk.Frame(self)

        container.pack(side="top", fill="both", expand=True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (StartPage, OligosPage):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="NSEW")

        self.show_frame(StartPage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()

class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        #self.master.title("PathoFinder Oligo DB")
        
        label = tk.Label(self, text="Start Page")
        label.grid(pady=10, padx=10)

        button1 = tk.Button(self, text="Oligos",
                         command=lambda:controller.show_frame(OligosPage))
        button1.grid()

class OligosPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        label = tk.Label(self, text="Start Page")
        label.grid(pady=10, padx=10)

        button1 = tk.Button(self, text="Back to Start Page",
                         command=lambda:controller.show_frame(StartPage))
        button1.grid()



if __name__ == "__main__":
    app = OligoDatabase()
    app.mainloop()
