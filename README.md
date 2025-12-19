Matthieu Beley - Class E design
V1.0 - 19/12/2025

Analytical equations describing the optimal ideal class E inverter are written in the form of a matlab code in "classE_design.m" for reproducibility.
Document "ClassE_design.pdf" explains its structure.

Charts are support for the design of such circuits. They allow for an exploration of the design space and help the designer to identify one or several design points.
Charts for the evolution of Class E key parameters are provided in the document "ClassE_Charts.pdf" as well.

Yet to be implemented --> analytical expressions for :
1) non-zero voltage slope at turn-on (non-ZdVS)
2) non-zero voltage drop at turn-on (non-ZVS)
3) Impact of stray inductance in series with the switch
4) Series filter tuned on harmonic n (class E frequency multiplier)

For information, the analytical expressions for the class E rectifier are the same, except for the reactance X which is of an opposite sign.
