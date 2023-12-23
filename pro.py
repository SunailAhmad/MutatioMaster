import tkinter as tk
from tkinter import ttk, filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib import pyplot as plt

class MutationAnalyzer:
    def __init__(self, reference_seq, mutated_seq):
        self.reference_seq = reference_seq
        self.mutated_seq = mutated_seq
        self.snps = self.snps_finder()
        self.insertions, self.deletions = self.insertion_deletion_finder()
        self.A, self.T, self.G, self.C = self.mutation_frequency_finder()

    def snps_finder(self):
        snps = []
        for i in range(len(self.reference_seq)):
            if self.reference_seq[i] != self.mutated_seq[i]:
                snps.append([i + 1, self.reference_seq[i], self.mutated_seq[i]])
        return snps

    def insertion_deletion_finder(self):
        insertions = []
        deletions = []
        if len(self.reference_seq) < len(self.mutated_seq):
            for i in range(len(self.reference_seq), len(self.mutated_seq)):
                insertions.append(self.mutated_seq[i])
        elif len(self.reference_seq) > len(self.mutated_seq):
            for i in range(len(self.mutated_seq), len(self.reference_seq)):
                deletions.append(self.reference_seq[i])
        return insertions, deletions

    def mutation_frequency_finder(self):
        A_counts, T_counts, G_counts, C_counts = 0, 0, 0, 0
        for i in self.snps:
            base_at_index_2 = i[2]
            if base_at_index_2 == 'A':
                A_counts += 1
            elif base_at_index_2 == 'T':
                T_counts += 1
            elif base_at_index_2 == 'G':
                G_counts += 1
            elif base_at_index_2 == 'C':
                C_counts += 1
            else:
                continue
        return A_counts, T_counts, G_counts, C_counts

    def visualize_mutation_frequency(self):
        # Plot the bar chart
        bases = ['A', 'T', 'G', 'C']
        counts = [self.A, self.T, self.G, self.C]

        plt.bar(bases, counts, color=['red', 'blue', 'green', 'purple'])
        plt.xlabel('Bases')
        plt.ylabel('Counts')
        plt.title('Mutation Frequency Histogram')
        plt.show()

        # Check if there are mutations before creating the pie chart
        if any(count > 0 for count in counts):
            # Create a pie chart
            sizes = counts
            colors = ['red', 'blue', 'green', 'purple']

            plt.pie(sizes, labels=bases, colors=colors, autopct='%1.1f%%', startangle=140)
            plt.axis('equal')  # Equal aspect ratio ensures that the pie is drawn as a circle.
            plt.title('Mutation Frequency Pie Chart')
            plt.show()
        else:
            print("No mutations to visualize.")

    def display_summary(self):
        print("Mutation Summary:")
        print(f"Number of SNPs: {len(self.snps)}")
        print(f"Insertions: {self.insertions}")
        print(f"Deletions: {self.deletions}")
        print(f"Mutation Frequencies: A={self.A}, T={self.T}, G={self.G}, C={self.C}")

    def highlight_mutations(self):
        highlighted_seq = list(self.reference_seq)
        for snp in self.snps:
            index = snp[0] - 1
            highlighted_seq[index] = f"({snp[1]}>{snp[2]})"
        return "".join(highlighted_seq)

    def save_results_to_file(self, filename="mutation_analysis_results.txt"):
        with open(filename, "w") as file:
            file.write("Mutation Analysis Results:\n")
            file.write("Reference Sequence: {}\n".format(self.reference_seq))
            file.write("Mutated Sequence: {}\n".format(self.mutated_seq))
            file.write("\nMutation Summary:\n")
            file.write("Number of SNPs: {}\n".format(len(self.snps)))
            file.write("Insertions: {}\n".format(self.insertions))
            file.write("Deletions: {}\n".format(self.deletions))
            file.write("Mutation Frequencies: A={}, T={}, G={}, C={}\n".format(self.A, self.T, self.G, self.C))
            file.write("\nMutated Sequence with Highlighted Mutations:\n")
            file.write("{}\n".format(self.highlight_mutations()))
        print("Results saved to {}".format(filename))

    def generate_detailed_report(self):
        print("Detailed Mutation Report:")
        for snp in self.snps:
            position, reference_base, mutated_base = snp
            print(f"Position: {position}, Reference Base: {reference_base}, Mutated Base: {mutated_base}")

    def visualize_mutation_positions(self):
        fig, ax = plt.subplots(figsize=(10, 1))
        ax.scatter([snp[0] for snp in self.snps], [0] * len(self.snps), marker='x', color='red')
        ax.set_yticks([])
        ax.set_xlabel('Sequence Position')
        ax.set_title('Mutation Positions on the Sequence')
        plt.show()

    def load_sequence_from_file(self, entry_var):
        file_path = filedialog.askopenfilename(title="Select a sequence file", filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if file_path:
            with open(file_path, "r") as file:
                sequence = file.read()
                entry_var.set(sequence)

    def summary_statistics(self):
        # Additional summary statistics
        num_mutations = len(self.snps)
        avg_distance_between_mutations = sum([self.snps[i][0] - self.snps[i - 1][0] for i in range(1, num_mutations)]) / (
                    num_mutations - 1) if num_mutations > 1 else 0
        mutation_types = [snp[1] + ">" + snp[2] for snp in self.snps]
        most_common_mutation_type = max(set(mutation_types), key=mutation_types.count)

        print("Summary Statistics:")
        print(f"Number of Mutations: {num_mutations}")
        print(f"Average Distance Between Mutations: {avg_distance_between_mutations:.2f} positions")
        print(f"Most Common Mutation Type: {most_common_mutation_type}")

class MutationAnalyzerGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Mutation Analyzer GUI")
        self.geometry("800x600")
        self.configure(bg='')  # Set background color

        self.reference_seq_var = tk.StringVar()
        self.mutated_seq_var = tk.StringVar()

        self.result_label_var = tk.StringVar()
        self.result_label_var.set("Results will be displayed here.")

        self.create_widgets()

    def create_widgets(self):
        # Entry widgets for user input
        reference_entry = ttk.Entry(self, textvariable=self.reference_seq_var, width=80, font=('Arial', 12))
        reference_entry.pack(pady=10, padx=10, side=tk.TOP, fill=tk.X)
        reference_entry.insert(0, "Enter reference sequence")

        mutated_entry = ttk.Entry(self, textvariable=self.mutated_seq_var, width=80, font=('Arial', 12))
        mutated_entry.pack(pady=10, padx=10, side=tk.TOP, fill=tk.X)
        mutated_entry.insert(0, "Enter mutated sequence")

        # Analyze Button
        analyze_button = ttk.Button(self, text="Analyze Sequences", command=self.analyze_sequences, style="TButton")
        analyze_button.pack(pady=10)

        # Load Sequences Button
        load_sequences_button = ttk.Button(self, text="Load Sequences from File", command=self.load_sequences_from_file, style="TButton")
        load_sequences_button.pack(pady=10)

        # Save Results Button
        save_results_button = ttk.Button(self, text="Save Results to File", command=self.save_results_to_file, style="TButton")
        save_results_button.pack(pady=10)

        # Additional Action Button
        additional_button = ttk.Button(self, text="Additional Action", command=self.additional_action, style="TButton")
        additional_button.pack(pady=10)

        # Matplotlib Figures for Mutation Frequency
        self.mutation_frequency_figure = Figure(figsize=(5, 4), dpi=100)
        self.mutation_frequency_canvas = FigureCanvasTkAgg(self.mutation_frequency_figure, master=self)
        self.mutation_frequency_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Matplotlib Figures for Mutation Positions
        self.mutation_positions_figure = Figure(figsize=(5, 4), dpi=100)
        self.mutation_positions_canvas = FigureCanvasTkAgg(self.mutation_positions_figure, master=self)
        self.mutation_positions_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Matplotlib Figures for Summary Statistics
        self.summary_statistics_figure = Figure(figsize=(5, 4), dpi=100)
        self.summary_statistics_canvas = FigureCanvasTkAgg(self.summary_statistics_figure, master=self)
        self.summary_statistics_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Navigation Toolbar for Matplotlib Figures
        self.toolbar = NavigationToolbar2Tk(self.mutation_frequency_canvas, self)
        self.toolbar.update()
        self.mutation_frequency_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Result Label
        result_label = tk.Label(self, textvariable=self.result_label_var, font=('Arial', 12), wraplength=700, justify='left')
        result_label.pack(pady=10)

    def analyze_sequences(self):
        reference_seq = self.reference_seq_var.get()
        mutated_seq = self.mutated_seq_var.get()

        if not reference_seq or not mutated_seq:
            self.result_label_var.set("Please enter both reference and mutated sequences.")
            return

        analyzer = MutationAnalyzer(reference_seq, mutated_seq)
        analyzer.display_summary()

        # Plot Mutation Frequency
        self.mutation_frequency_figure.clear()
        analyzer.visualize_mutation_frequency()
        self.mutation_frequency_canvas.draw()

        # Plot Mutation Positions
        self.mutation_positions_figure.clear()
        analyzer.visualize_mutation_positions()
        self.mutation_positions_canvas.draw()

        # Plot Summary Statistics
        self.summary_statistics_figure.clear()
        analyzer.summary_statistics()
        self.summary_statistics_canvas.draw()

        # Update result label
        result_text = f"Mutation Summary:\nNumber of SNPs: {len(analyzer.snps)}\nInsertions: {analyzer.insertions}\n" \
                      f"Deletions: {analyzer.deletions}\nMutation Frequencies: A={analyzer.A}, T={analyzer.T}, " \
                      f"G={analyzer.G}, C={analyzer.C}"
        self.result_label_var.set(result_text)

        return analyzer

    def load_sequences_from_file(self):
        # Load reference sequence from file
        self.load_sequence_from_file(self.reference_seq_var)

        # Load mutated sequence from file
        self.load_sequence_from_file(self.mutated_seq_var)

    def save_results_to_file(self):
        analyzer = self.analyze_sequences()  # Create an analyzer instance to access results
        save_path = filedialog.asksaveasfilename(title="Save Results", defaultextension=".txt",
                                                  filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if save_path:
            analyzer.save_results_to_file(filename=save_path)

    def load_sequence_from_file(self, entry_var):
        file_path = filedialog.askopenfilename(title="Select a sequence file",
                                                filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if file_path:
            with open(file_path, "r") as file:
                sequence = file.read()
                entry_var.set(sequence)

    def additional_action(self):
        # Add your custom functionality here
        print("Performing additional action")

# Create and run the application
if __name__ == "__main__":
    app = MutationAnalyzerGUI()
    app.mainloop()
