# Import libraries
import time
import tkinter as tk
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk
import subprocess
import os

# Global variables
new_bacteria = None
database_path = None
selected_treshold = None
loading_bar = None
download_final = None
download_all = None
link_label = None
result_widget = None

## Create root
root = tk.Tk()
root.title("BACTERIAL IDENTIFICATION")

def button1_click():
    global new_bacteria
    new_bacteria = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta")])
    if new_bacteria:
        print(new_bacteria)
        lbl5 = Label(root, text="Successfully upload", fg='white', font=("Helvetica", 14), bg="cornflowerblue")
        lbl5.place(x=520, y=55)

def button2_click():
    global database_path
    database_path = filedialog.askdirectory()
    if database_path:
        print(database_path)
        lbl6 = Label(root, text="Successfully upload", fg='white', font=("Helvetica", 14), bg="cornflowerblue")
        lbl6.place(x=520, y=85)

def button3_click():
    global new_bacteria, database_path, selected_threshold, download_final, download_all, link_label
    if new_bacteria and database_path and (cd_hit.get()==1 or fastAni.get()==1 or blast.get()==1):
        button3.config(text="Runing analysis...", bg="green")
        root.update_idletasks()

        if download_final:
            download_final.destroy()
        if download_all:
            download_all.destroy()
        if link_label:
            link_label.place_forget()
        if result_widget:
            result_widget.place_forget()

        # Loading bar
        loading_bar = ttk.Progressbar(root, mode='indeterminate')
        loading_bar.place(x=350, y=165, width=300)

        def animate_loading_bar():
            for i in range(100):
                loading_bar.step(1)
                root.update_idletasks()
                root.after(10)

        # Starting loading bar animation
        animate_loading_bar()

        if cd_hit.get()==1 and blast.get()==0:
            methods = "cd_hit"
        elif cd_hit.get()==0 and blast.get()==1:
            methods = "blast"
        else:
            methods = "all"
        print(methods)

        ## Starting analysis (shell script)
        start_time = time.time()
        cmd = ["./pipeline.sh", new_bacteria, database_path, methods, str(selected_threshold)]
        subprocess.run(cmd)
        end_time = time.time()
        elapsed_time = end_time-start_time
        print("Doba analýzy: ", elapsed_time, " sekund")

        ## Stop analysis
        button3.config(text="Analysis finished!", bg="green")
        root.update_idletasks()

        # Stop loading_bar
        loading_bar.stop()
        loading_bar.destroy()
        analyze_results()
        show_download_links()

    # Errors
    elif not new_bacteria and not database_path:
        messagebox.showerror("Error", "Please load bacterial genome and reference database!")
    elif not new_bacteria:
        messagebox.showerror("Error", "You don't load unknown bacterial genome!")
    elif not database_path:
        messagebox.showerror("Error", "You don't load reference database!")
    elif cd_hit.get() == 0 and blast.get()== 0:
        messagebox.showerror("Error", "Select method!")


def analyze_results():
    global result_widget
    final_file = 'final.txt'
    all_file = "all_possible_bacteria.txt"

    all_file_exists = os.path.exists(final_file) and os.path.getsize(final_file) > 0
    final_file_exists = os.path.exists(all_file) and os.path.getsize(all_file) > 0

    if not final_file_exists and not all_file_exists:
        result_widget = Label(root, text="It's probably a completely new bacterium!", fg='white', font=("Helvetica", 16), bg="cornflowerblue")
        result_widget.place(x=100, y=230)
    else:
        result_widget = Label(root, text="Bacteria with significant similarity were found!", fg='white', font=("Helvetica", 16), bg="cornflowerblue")
        result_widget.place(x=100, y=220)


def show_download_links():
    global download_final, download_all, link_label
    link_label = Label(root, text="Download Result Files:", fg='black', font=("Helvetica", 12), bg="cornflowerblue")
    link_label.place(x=250, y=250)

    # Download links
    download_final = Button(root, text="Download final results", command=download_final_results)
    download_final.place(x=250, y=280)

    download_all = Button(root, text="Download all possible results", command=download_all_results)
    download_all.place(x=250, y=320)

def download_final_results():
    final_results_file = "final.txt"
    if os.path.exists(final_results_file):
        os.system(f"xdg-open {final_results_file}")
    else:
        messagebox.showerror("Error", "Final results file not found!")

def download_all_results():
    final_results_file = "all_possible_bacteria.txt"
    if os.path.exists(final_results_file):
        os.system(f"xdg-open {final_results_file}")
    else:
        messagebox.showerror("Error", "All possible results file not found!")


# Create IntVar for checkboxes
cd_hit = tk.IntVar()
fastAni = tk.IntVar()
blast = tk.IntVar()

# Create checkboxes
checkbox1 = tk.Checkbutton(root, text="Cd-hit", variable=cd_hit, bg="cornflowerblue")
checkbox2 = tk.Checkbutton(root, text="BLAST", variable=blast, bg="cornflowerblue")

# Create buttons
button1 = tk.Button(root, text="Upload", command=button1_click)
button2 = tk.Button(root, text="Upload", command=button2_click)
button3 = tk.Button(root, text="Start analysis!", command=button3_click)

# Create MenuButton
options = [95, 98, 99]
selected = IntVar()

menu_button = tk.Menubutton(root, text="Set threshold", indicatoron=True, borderwidth=1, relief="raised")
menu_button.pack()

menu = tk.Menu(menu_button, tearoff=False)
menu_button["menu"] = menu

for option in options:
    menu.add_radiobutton(label=option, variable=selected, value=option)

selected_threshold = selected.get()

# Set default threshold
if selected_threshold == 0:
    selected_threshold = 95

# Create label
lbl1 = Label(root, text="Please set parameters of the tool!", fg='red', font=("Helvetica", 20), bg="cornflowerblue")
lbl2 = Label(root, text="Select methods", fg='black', font=("Helvetica", 14), bg="cornflowerblue")
lbl3 = Label(root, text="Upload genome of the bacteria", fg='black', font=("Helvetica", 14), bg="cornflowerblue")
lbl4 = Label(root, text="Upload path to the database ", fg='black', font=("Helvetica", 14), bg="cornflowerblue")
lbl5 = Label(root, text="Set threshold", fg='black', font=("Helvetica", 14), bg="cornflowerblue")


# Wudget location in root
checkbox1.place(x=30, y=60)
checkbox2.place(x=30, y=80)
button1.place(x=440, y=55)
button2.place(x=440, y=85)
button3.place(x=200, y=160)
lbl1.place(x=0, y=0)
lbl2.place(x=0, y=40)
lbl3.place(x=180, y=60)
lbl4.place(x=180, y=90)
lbl5.place(x=0, y=130)
menu_button.place(x=30, y=160)

# Size of widget
root.geometry("700x400")
root.config(bg="cornflowerblue")

root.mainloop()
