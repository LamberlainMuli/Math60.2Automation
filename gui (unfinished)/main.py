import tkinter as tk
from tkinter import messagebox
import math
from scipy.stats import norm
from stat_utility import *

# GUI code starts here
class MainMenu(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Main Menu")
        tk.Button(self, text="Sampling Distribution of the Sample Mean", command=self.open_mean_page).pack()
        tk.Button(self, text="Sampling Distribution of the Sample Proportion", command=self.open_proportion_page).pack()

    def open_mean_page(self):
        MeanPage()

    def open_proportion_page(self):
        ProportionPage()


class MeanPage(tk.Toplevel):
    def __init__(self):
        super().__init__()
        self.title("Sampling Distribution of the Sample Mean")
        self.mean = tk.StringVar()
        self.sd = tk.StringVar()
        self.n = tk.StringVar()
        self.lower_bound = tk.StringVar()
        self.upper_bound = tk.StringVar()

        tk.Label(self, text="Mean").grid(row=0, column=0)
        tk.Entry(self, textvariable=self.mean).grid(row=0, column=1)
        tk.Label(self, text="Standard Deviation").grid(row=1, column=0)
        tk.Entry(self, textvariable=self.sd).grid(row=1, column=1)
        tk.Label(self, text="Sample Size").grid(row=2, column=0)
        tk.Entry(self, textvariable=self.n).grid(row=2, column=1)
        tk.Label(self, text="Lower Bound").grid(row=3, column=0)
        tk.Entry(self, textvariable=self.lower_bound).grid(row=3, column=1)
        tk.Label(self, text="Upper Bound").grid(row=4, column=0)
        tk.Entry(self, textvariable=self.upper_bound).grid(row=4, column=1)
        tk.Button(self, text="Calculate", command=self.calculate).grid(row=5, column=0, columnspan=2)

    def calculate(self):
        try:
            mean = float(self.mean.get())
            sd = float(self.sd.get())
            n = int(self.n.get())
            lower_bound = self.lower_bound.get()
            upper_bound = self.upper_bound.get()

            if lower_bound and upper_bound:
                prob = get_probability_between_mean(float(lower_bound), float(upper_bound), mean, sd, n)
            elif lower_bound:
                prob = get_probability_greater_than_mean(float(lower_bound), mean, sd, n)
            elif upper_bound:
                prob = get_probability_less_than_mean(float(upper_bound), mean, sd, n)
            else:
                prob = None

            sampling_distribution = f"N({mean}, {get_standard_error_mean(sd, n)})"
            messagebox.showinfo("Result", f"Sampling Distribution: {sampling_distribution}\nProbability: {prob}")
        except ValueError:
            messagebox.showerror("Error", "Invalid input. Please enter numbers only.")


class ProportionPage(tk.Toplevel):
    def __init__(self):
        super().__init__()
        self.title("Sampling Distribution of the Sample Proportion")
        self.p = tk.StringVar()
        self.n = tk.StringVar()
        self.lower_bound = tk.StringVar()
        self.upper_bound = tk.StringVar()

        tk.Label(self, text="Proportion").grid(row=0, column=0)
        tk.Entry(self, textvariable=self.p).grid(row=0, column=1)
        tk.Label(self, text="Sample Size").grid(row=1, column=0)
        tk.Entry(self, textvariable=self.n).grid(row=1, column=1) 
        tk.Label(self, text="Lower Bound").grid(row=2, column=0)
        tk.Entry(self, textvariable=self.lower_bound).grid(row=2, column=1)
        tk.Label(self, text="Upper Bound").grid(row=3, column=0)
        tk.Entry(self, textvariable=self.upper_bound).grid(row=3, column=1)
        tk.Button(self, text="Calculate", command=self.calculate).grid(row=4, column=0, columnspan=2)

    def calculate(self):
        try:
            p = float(self.p.get())
            n = int(self.n.get())
            lower_bound = self.lower_bound.get()
            upper_bound = self.upper_bound.get()

            if lower_bound and upper_bound:
                prob = get_probability_between_proportion(float(lower_bound), float(upper_bound), p, n)
            elif lower_bound:
                prob = get_probability_greater_than_proportion(float(lower_bound), p, n)
            elif upper_bound:
                prob = get_probability_less_than_proportion(float(upper_bound), p, n)
            else:
                prob = None

            sampling_distribution = f"N({p}, {get_standard_error_proportion(n, p)})"
            messagebox.showinfo("Result", f"Sampling Distribution: {sampling_distribution}\nProbability: {prob}")
        except ValueError:
            messagebox.showerror("Error", "Invalid input. Please enter numbers only.")


if __name__ == "__main__":
    MainMenu().mainloop()