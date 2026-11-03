Matthieu Beley - Class E design

Updated - 03/11/2026

#
Analytical equations describing the optimal ideal class E inverter are written in the form of a matlab code for reproducibility.
The script "Generic_Class_E_design.m" implements the analytical expressions with variable parameters :
* q --> the input LC tank resonant frequency ratio
* D --> the ON duty cycle
* K --> the switch voltage at turn-on (0 if ZVS)
* s --> the switch voltage derivative at turn-on (0 is ZdVS, or ZCS)

Document "ClassE_design.pdf" explains its structure.

#
Charts are support for the design of such circuits. They allow for an exploration of the design space and help the designer to identify one or several design points.
Charts for the evolution of Class E key parameters are provided in the document "ClassE_Charts.pdf" as well.

#
Yet to be implemented --> analytical expressions for :
1) Impact of stray inductance in series with the switch
2) Series filter tuned on harmonic n (class E frequency multiplier)

#
For information, the analytical expressions for the class E rectifier are the same, except for the reactance X which is of an opposite sign.
