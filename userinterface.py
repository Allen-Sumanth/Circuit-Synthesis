import tkinter as tk
import computation as comp

# Function to handle the "synthesize" button click
def synthesize():
    numerator = numerator_entry.get()
    denominator = denominator_entry.get()
    synthesis_format = synthesis_format_var.get()
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")
    comp.synthesize(numerator, denominator)
    
    root.destroy()

options = ["Foster 1", "Foster 2", "Cauer 1", "Cauer 2"]

def update_option_index(*args):
    comp.optno = options.index(synthesis_format_var.get())
    print(f"Option index: {comp.optno}")
    
root = tk.Tk()
root.title("Network Synthesis Input")

# Label and input field for the numerator
numerator_label = tk.Label(root, text="Enter the numerator of the Laplace function:\nExample: 2*s**2 + 3*s + 6")
numerator_label.pack(pady=5, padx=15)
numerator_entry = tk.Entry(root, width=50)
numerator_entry.pack(pady=5, padx=15)

# Label and input field for the denominator
denominator_label = tk.Label(root, text="Enter the denominator of the Laplace function:\nExample: 2*s**2 + 3*s + 6")
denominator_label.pack(pady=5)
denominator_entry = tk.Entry(root, width=50)
denominator_entry.pack(pady=5)

# Label and dropdown menu for the synthesis format
synthesis_format_label = tk.Label(root, text="Synthesis Format")
synthesis_format_label.pack(pady=5)
synthesis_format_var = tk.StringVar(root)
synthesis_format_var.set("Foster 1")  # default value
synthesis_format_var.trace("w", update_option_index)
synthesis_format_menu = tk.OptionMenu(root, synthesis_format_var, *options)
synthesis_format_menu.pack(pady=5)

synthesize_button = tk.Button(root, text="Synthesize", command=synthesize)
synthesize_button.pack(pady=20)
