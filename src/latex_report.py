import os
import subprocess

class LatexReport:
    def __init__(self, filename="report.tex", image_path="image.png"):
        self.filename = filename
        self.image_path = image_path
        self.latex_content = ""

    def write_latex(self):
        """
        Write the LaTeX content for a simple report with an image.
        """
        self.latex_content = rf"""
        \documentclass{{article}}
        \usepackage{{amsmath}}
        \usepackage{{graphicx}}
        \title{{Sample LaTeX Report with Image}}
        \author{{Your Name}}
        \date{{\today}}

        \begin{{document}}

        \maketitle

        \section{{Introduction}}
        This is an example of a simple LaTeX report generated using Python.
        
        \section{{Mathematics}}
        Here is a simple equation:
        \begin{{equation}}
            E = mc^2
        \end{{equation}}

        \section{{Image Section}}
        The image below is included in this report:

        \begin{{center}}
        \includegraphics[width=1\textwidth]{{{self.image_path}}}
        \end{{center}}

        \section{{Conclusion}}
        We have successfully generated a LaTeX report with an image using Python and compiled it into a PDF.

        \end{{document}}
        """
        # Write the LaTeX content to a .tex file
        with open(self.filename, 'w') as f:
            f.write(self.latex_content)
        print(f"LaTeX file '{self.filename}' written.")

    def compile_latex(self):
        """
        Compile the LaTeX file to generate a PDF using pdflatex.
        """
        try:
            subprocess.run(['pdflatex', self.filename], check=True)
            print("LaTeX compiled successfully!")
        except subprocess.CalledProcessError as e:
            print(f"Error during LaTeX compilation: {e}")
        finally:
            # Clean up auxiliary files generated during compilation
            for ext in ['.aux', '.log', '.out']:
                aux_file = self.filename.replace('.tex', ext)
                if os.path.exists(aux_file):
                    os.remove(aux_file)

if __name__ == "__main__":
    # Ensure you have the path to your image file
    image_path = r"C:/Users/Jason/Desktop/stochICE/plot.pdf"  # Replace with the path to your image file

    # Create an instance of the LatexReport class
    report = LatexReport(image_path=image_path)

    # Write the LaTeX content to a file
    report.write_latex()

    # Compile the LaTeX file to produce a PDF
    report.compile_latex()
