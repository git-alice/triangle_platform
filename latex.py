import sympy as sym
import sympy.printing as printing
# delta__y_l = sym.symbols('Delta__y_l')
# print(printing.latex(delta__y_l))

def print_tex(data):
    return '''
\\title{A Very Simple \\LaTeXe{} Template}
\\author{
        Vitaly Surazhsky \\\\
                Department of Computer Science\\\\
        Technion---Israel Institute of Technology\\\\
        Technion City, Haifa 32000, \\underline{Israel}
            \\and
        Yossi Gil\\\\
        Department of Computer Science\\\\
        Technion---Israel Institute of Technology\\\\
        Technion City, Haifa 32000, \\underline{Israel}
}
\\date{\\today}

\\documentclass[12pt]{article}

\\begin{document}
$$%s$$
\\end{document}
''' % str(printing.latex(data))