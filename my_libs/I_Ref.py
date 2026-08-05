import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci

from my_libs import Constants as const
from my_libs import T_E as tecore
from my_libs import Probs as prob
from my_libs import Quantum_Core as qcore
from my_libs import TMM as tmm

## Constants ##
a = const.a ## Fine structure constant
c = const.c ## Speed of light in a vacuum
hb = const.hb ## reduced plancks constant
h = const.h ## Plancks constant
me = const.me ##electron mass (kg)
q = const.q ## electron charge
kb = const.kb ## Boltzmans constant

## Integral limits ##
int_low = const.E_LOW
int_high = const.E_HIGH

def Cal_TE_Alumina():
    T = 300
    Area_1 = 1e-4
    
    params_45 = {
        "EA_Cond": 4.45*q,
        "EA_Barrier": const.EA_Al2O3,
        "m*_Barrier": 0.2,
        "Bar_Thickness": 4.5e-9
    }
    params_50 = {
        "EA_Cond": 4.45*q,
        "EA_Barrier": const.EA_Al2O3,
        "m*_Barrier": 0.2,
        "Bar_Thickness": 5.0e-9
    }
    params_54 = {
        "EA_Cond": 4.45*q,
        "EA_Barrier": const.EA_Al2O3,
        "m*_Barrier": 0.2,
        "Bar_Thickness": 5.4e-9
    }
    params_59 = {
        "EA_Cond": 4.45*q,
        "EA_Barrier": const.EA_Al2O3,
        "m*_Barrier": 0.2,
        "Bar_Thickness": 5.9e-9
    }
    
    Voltage_1 = np.linspace(0.01, 10, 10000)

    J45 = lambda E, V: tecore.N(E, (V)*q, 0, T)*prob.Barrier_Prob(V, E, params_45)
    J50 = lambda E, V: tecore.N(E, (V)*q, 0, T)*prob.Barrier_Prob(V, E, params_50)
    J54 = lambda E, V: tecore.N(E, (V)*q, 0, T)*prob.Barrier_Prob(V, E, params_54)
    J59 = lambda E, V: tecore.N(E, (V)*q, 0, T)*prob.Barrier_Prob(V, E, params_59)
    
    I_45 = []
    I_50 = []
    I_54 = []
    I_59 = []
    
    for val in Voltage_1:
        result, error = sci.quad(J45, (val*0.45 - int_low)*q, (val*0.45 + int_high)*q, args = (val*0.45,), epsabs=1e-12, epsrel=1e-12)
        I_45.append(result*tecore.Constant(0.2)*Area_1)
        result, error = sci.quad(J50, (val*0.5 - int_low)*q, (val*0.5 + int_high)*q, args = (val*0.5,), epsabs=1e-12, epsrel=1e-12)
        I_50.append(result*tecore.Constant(0.2)*Area_1)
        result, error = sci.quad(J54, (val*0.54 - int_low)*q, (val*0.54 + int_high)*q, args = (val*0.54,), epsabs=1e-12, epsrel=1e-12)
        I_54.append(result*tecore.Constant(0.2)*Area_1)
        result, error = sci.quad(J59, (val*0.59 - int_low)*q, (val*0.59 + int_high)*q, args = (val*0.59,), epsabs=1e-12, epsrel=1e-12)
        I_59.append(result*tecore.Constant(0.2)*Area_1)
        
    x45 = [0.021231422505307854, 0.7643312101910827, 1.592356687898089, 2.505307855626327, 3.906581740976645, 5.095541401273885, 6.645435244161359, 7.558386411889596, 8.768577494692144, 9.872611464968152]
    y45 = [2.1564571420575197e-10, 1.2490419281239316e-9, 1.2992476247809846e-8, 1.9486374793329706e-7, 0.000002026963674778543, 0.000016928027542502107, 0.00009804885869815224, 0.0004559573337078687, 0.003807891773914453, 0.049334875234468986]
    
    x50 = [0.021231422505307854, 1.1464968152866242, 2.3991507430997876, 3.5881104033970272, 4.861995753715498, 6.199575371549893, 6.985138004246284, 7.579617834394904, 8.51380042462845, 9.341825902335456, 9.872611464968152]
    y50 = [1.3900543618179395e-10, 1.002818861742005e-9, 7.78391426873118e-9, 8.096791546732073e-8, 6.76197171239162e-7, 0.000005444290873268693, 0.000037865143653810273, 0.00031622776601683794, 0.002281343644384445, 0.014217123497267647, 0.04585315181043681]
    
    x54 = [-0.021231422505307854, 0.9978768577494691, 2.314225053078556, 3.4394904458598723, 4.501061571125265, 5.329087048832271, 6.13588110403397, 6.7515923566878975, 7.346072186836517, 8.280254777070063, 8.9171974522293, 9.745222929936306]
    y54 = [1.2919536813314063e-10, 8.985561733646273e-10, 2.793896212615348e-9, 1.2992476247809846e-8, 7.805856046866354e-8, 3.499530475916605e-7, 0.0000027163475124538404, 0.000019596388424494337, 0.00014664226793541493, 0.001057912814445164, 0.006355920905444216, 0.03818625694353637]
    
    x59 = [0.25477707006369427, 2.08067940552017, 3.2908704883227173, 4.522292993630573, 5.350318471337579, 5.796178343949045, 6.518046709129512, 7.133757961783439, 8.046709129511676, 8.938428874734607, 9.681528662420382]
    y59 = [1.2919536813314063e-10, 8.051336362871474e-10, 2.009918668908342e-9, 1.2075554725799471e-8, 9.036291245156562e-8, 5.233911873236153e-7, 0.000007295930437576438, 0.00006800147569615066, 0.0009479220291588477, 0.005102981112326089, 0.04756216483377386]
    
    plt.plot(Voltage_1, np.log10(I_59), label = "5.9nm Calc")
    plt.plot(Voltage_1, np.log10(I_54), label = "5.4nm Calc")
    plt.plot(Voltage_1, np.log10(I_50), label = "5.0nm Calc")
    plt.plot(Voltage_1, np.log10(I_45), label = "4.5nm Calc")
    
    plt.plot(x59, np.log10(y59), label = '5.9nm exp', linestyle = '--')
    plt.plot(x54, np.log10(y54), label = '5.4nm exp', linestyle = '--')
    plt.plot(x50, np.log10(y50), label = '5.0nm exp', linestyle = '--')
    plt.plot(x45, np.log10(y45), label = '4.5nm exp', linestyle = '--')
    
    plt.legend()
    plt.title("Tsu-Esaki Confirmation (TiN/Al2O3/TiN)")
    plt.xlabel("E-field (MV/cm)")
    plt.ylabel("Current A/cm^2")
    plt.ylim([-10, 1])
    plt.show()

def Cal_TE_BN():
    T = 300
    Area_2 = (1e-6)**2 * 1e6
    
    params_1 = {
        "EA_Cond": 4.6*q,
        "EA_Barrier": const.EA_hBN_TYP,
        "m*_Barrier": 0.5,
        "Bar_Thickness": 1*const.ML_hBN
    }
    params_2 = {
        "EA_Cond": 4.6*q,
        "EA_Barrier": const.EA_hBN_TYP,
        "m*_Barrier": 0.5,
        "Bar_Thickness": 2*const.ML_hBN
    }
    params_3 = {
        "EA_Cond": 4.6*q,
        "EA_Barrier": const.EA_hBN_TYP,
        "m*_Barrier": 0.5,
        "Bar_Thickness": 3*const.ML_hBN
    }
    params_4 = {
        "EA_Cond": 4.6*q,
        "EA_Barrier": const.EA_hBN_TYP,
        "m*_Barrier": 0.5,
        "Bar_Thickness": 4*const.ML_hBN
    }
    
    Voltage_2 = np.linspace(0.01, 1, 10000)
    
    J_1 = lambda E, V: tecore.N(E, (V)*q, 0, T)*prob.Barrier_Prob(V, E, params_1)
    J_2 = lambda E, V: tecore.N(E, (V)*q, 0, T)*prob.Barrier_Prob(V, E, params_2)
    J_3 = lambda E, V: tecore.N(E, (V)*q, 0, T)*prob.Barrier_Prob(V, E, params_3)
    J_4 = lambda E, V: tecore.N(E, (V)*q, 0, T)*prob.Barrier_Prob(V, E, params_4)
    
    I_1 = []
    I_2 = []
    I_3 = []
    I_4 = []
    
    for val in Voltage_2:
        result, error = sci.quad(J_1, (val - int_low)*q, (val + int_high)*q, args = (val,), epsabs=1e-12, epsrel=1e-12)
        I_1.append(result*tecore.Constant(0.5)*Area_2)
        result, error = sci.quad(J_2, (val - int_low)*q, (val + int_high)*q, args = (val,), epsabs=1e-12, epsrel=1e-12)
        I_2.append(result*tecore.Constant(0.5)*Area_2)
        result, error = sci.quad(J_3, (val - int_low)*q, (val + int_high)*q, args = (val,), epsabs=1e-12, epsrel=1e-12)
        I_3.append(result*tecore.Constant(0.5)*Area_2)
        result, error = sci.quad(J_4, (val - int_low)*q, (val + int_high)*q, args = (val,), epsabs=1e-12, epsrel=1e-12)
        I_4.append(result*tecore.Constant(0.5)*Area_2)
        
    x4 = [0.014925373134328358, 0.03781094527363184, 0.11243781094527362, 0.1781094527363184, 0.2825870646766169, 0.46467661691542284, 0.582089552238806, 0.7044776119402985, 0.8079601990049751, 0.9393034825870646, 0.9920398009950249]
    y4 = [0.0009378028619870069, 0.004889258820094655, 0.009464455149705116, 0.015819880086366582, 0.023686559963647765, 0.042607109440588355, 0.07121797837436074, 0.10279031364666079, 0.1404139829332773, 0.1989778715200249, 0.210237006116749]
    
    x3 = [0.0059701492537313425, 0.03980099502487562, 0.1263681592039801, 0.21592039800995025, 0.32238805970149254, 0.4338308457711442, 0.5592039800995025, 0.6686567164179104, 0.7810945273631841, 0.9064676616915422, 0.9880597014925373]
    y3 = [0.005072000443730913, 0.05714419400520047, 0.1989778715200249, 0.4460693572977468, 0.7456069716673378, 1.269362480415192, 1.900572592447801, 2.793927545571262, 3.482031289161529, 5.025684476172322, 5.820231231507102]
    
    x2 = [0.003980099502487562, 0.051741293532338306, 0.16517412935323383, 0.26865671641791045, 0.36915422885572136, 0.481592039800995, 0.5781094527363184, 0.6965174129353233, 0.8069651741293532, 0.9512437810945273, 0.9960199004975123]
    y2 = [0.6093412200099704, 3.418727944052855, 15.11064033315543, 28.19675970124227, 47.131012851980756, 73.20518334381646, 91.23455593400381, 120.13824264608023, 167.15045753964486, 212.1745193850578, 224.1803844267449]
    
    x1 = [-0.0009950248756218905, 0.031840796019900496, 0.10149253731343283, 0.20597014925373133, 0.3243781094527363, 0.4497512437810945, 0.56318407960199, 0.6805970149253731, 0.772139303482587, 0.8736318407960199, 0.9482587064676616, 0.9950248756218905]
    y1 = [4.934317508910219, 50.72000443730913, 161.12811066858086, 410.71996728355685, 855.6002766702244, 1483.5928821285838, 2026.6227258262868, 2668.669675594201, 3325.9242264836835, 3923.064947006565, 4800.372071897567, 4979.791451962157]
    
    plt.plot(Voltage_2, np.log10(I_4), label = "1.38nm (4 Layers)")
    plt.plot(Voltage_2, np.log10(I_3), label = "1.02nm (3 Layers)")
    plt.plot(Voltage_2, np.log10(I_2), label = "0.68nm (2 Layers)")
    plt.plot(Voltage_2, np.log10(I_1), label = "0.34nm (1 Layers)")
    
    plt.plot(x4, np.log10(y4), label = '4 layers exp', linestyle = '--')
    plt.plot(x3, np.log10(y3), label = '3 layers exp', linestyle = '--')
    plt.plot(x2, np.log10(y2), label = '2 layers exp', linestyle = '--')
    plt.plot(x1, np.log10(y1), label = '1 layers exp', linestyle = '--')
    
    plt.legend()
    plt.title("Tsu-Esaki Confirmation (Gr/h-BN/Gr)")
    plt.xlabel("Volts")
    plt.ylabel("Current uA/um^2")
    plt.show()