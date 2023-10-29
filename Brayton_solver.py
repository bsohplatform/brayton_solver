from dataclasses import dataclass, field
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

@dataclass
class Fluid_stream:
    fluid:str
    T:float = 0 
    p:float = 0 
    h:float = 0 
    m:float = 0
    q:float = 0
    s:float = 0
    
    
@dataclass
class Settings:
    hot_dp:float = 0.0
    hot_T_pp:float = 0.0 
    cold_dp:float = 0.0
    cold_T_pp:float = 0.0
    
    ihx_cold_dp:float = 0.0
    ihx_hot_dp:float = 0.0
    ihx_eff:float = 0.0
    
    comp_eff:float = 0.0
    expand_eff:float = 0.0
    mech_eff:float = 0.95
    gen_eff: float = 0.99
    
    N_element:int = 30
    
@dataclass
class Outputs:
    COP:float = 0.0
    comp_W:float = 0.0
    expand_W:float = 0.0
    ihx_Q:float = 0.0
    hot_Q:float = 0.0
    cold_Q:float = 0.0
    hot_eff:float = 0.0
    cold_eff:float = 0.0
    T_cold_sec: list = field(default_factory=list)
    T_hot_sec: list = field(default_factory=list)

@staticmethod
class Aux:
    def TS_fluid_prop(x, fluid, index):
        if index == 'H':
            if fluid == 'na':
                y = 2386*x
            elif 'flibe':
                y = flibe_TH(x)
            elif 'solarsalt':
                y = 0.086*x**2+1396.0182*x
            else:
                print('Fluid is not defined')
        elif index == 'C':
            if fluid == 'na':
                y = 2386
            elif 'flibe':
                y = (1.6582-4.2395*0.0002*x+1.4847*0.0000003*x**2-2992.6*x**-2)*1e3
            elif 'solarsalt':
                y = 0.172*x+1396.0182
            else:
                print('Fluid is not defined')
        elif index == 'T':
            if fluid == 'na':
                y = x/2386
            elif 'flibe':
                y = fsolve(flibe_TH, 371)
            elif 'solarsalt':
                y = (-1396.0182+(1.948866814731240e+06+0.3440*x)**0.5)/0.1720;
            else:
                print('Fluid is not defined')
        
        def flibe_TH(x):
            y = (-365.77+1.6582*x-4.2395*0.0001*x**2+1.4847*0.0000001*x**3+2992.6*x**-1)*1e3
            
            return y
        
        return y

class Brayton():
    def __init__(self, tol=1.0e-4):
        self.tol = tol
        
    def Input_Processing(self, second):
        if second['hot_in'].p <= 0.0:
            second['hot_in'].p = 101300.0
        if second['hot_out'].p <= 0.0:
            second['hot_out'].p = 101300.0
        if second['hot_in'].T <= 0.0:
            process_case = 0
        if second['hot_out'].T <= 0.0:
            process_case = 1
        if second['hot_in'].m <= 0.0 and second['hot_out'].m <= 0.0:
            process_case = 2
        else:
            if second['hot_in'].m == 0:
                second['hot_in'].m = second['hot_out'].m
            else:
                second['hot_out'].m = second['hot_in'].m
        
        if second['cold_in'].p <= 0.0:
            second['cold_in'].p = 101300.0
        if second['cold_out'].p <= 0.0:
            second['cold_out'].p = 101300.0
        if second['cold_in'].T <= 0.0:
            process_case = 3
        if second['cold_out'].T <= 0.0:
            process_case = 4
        if second['cold_in'].m <= 0.0 and second['cold_out'].m <= 0.0:
            process_case = 5
        else:
            if second['cold_in'].m == 0:
                second['cold_in'].m = second['cold_out'].m
            else:
                second['cold_out'].m = second['cold_in'].m
                
        if process_case == 0:
            second['hot_out'].h = CP.PropsSI('H','T',second['hot_out'].T,'P',second['hot_out'].p,second['hot_out'].fluid)
            second['cold_in'].h = CP.PropsSI('H','T',second['cold_in'].T,'P',second['cold_in'].p,second['cold_in'].fluid)
            second['cold_out'].h = CP.PropsSI('H','T',second['cold_out'].T,'P',second['cold_out'].p,second['cold_out'].fluid)
            second['cold_in'].q = (second['cold_out'].h - second['cold_in'].h)*second['cold_in'].m
            second['cold_out'].q = second['cold_in'].q
        elif process_case == 1:
            second['hot_in'].h = CP.PropsSI('H','T',second['hot_in'].T,'P',second['hot_in'].p,second['hot_in'].fluid)
            second['cold_in'].h = CP.PropsSI('H','T',second['cold_in'].T,'P',second['cold_in'].p,second['cold_in'].fluid)
            second['cold_out'].h = CP.PropsSI('H','T',second['cold_out'].T,'P',second['cold_out'].p,second['cold_out'].fluid)
            second['cold_in'].q = (second['cold_out'].h - second['cold_in'].h)*second['cold_in'].m
            second['cold_out'].q = second['cold_in'].q
        elif process_case == 2:
            second['hot_in'].h = CP.PropsSI('H','T',second['hot_in'].T,'P',second['hot_in'].p,second['hot_in'].fluid)
            second['hot_out'].h = CP.PropsSI('H','T',second['hot_out'].T,'P',second['hot_out'].p,second['hot_out'].fluid)
            second['cold_in'].h = CP.PropsSI('H','T',second['cold_in'].T,'P',second['cold_in'].p,second['cold_in'].fluid)
            second['cold_out'].h = CP.PropsSI('H','T',second['cold_out'].T,'P',second['cold_out'].p,second['cold_out'].fluid)
            second['cold_in'].q = (second['cold_out'].h - second['cold_in'].h)*second['cold_in'].m
            second['cold_out'].q = second['cold_in'].q
        elif process_case == 3:
            second['hot_in'].h = CP.PropsSI('H','T',second['hot_in'].T,'P',second['hot_in'].p,second['hot_in'].fluid)
            second['hot_out'].h = CP.PropsSI('H','T',second['hot_out'].T,'P',second['hot_out'].p,second['hot_out'].fluid)
            second['cold_out'].h = CP.PropsSI('H','T',second['cold_out'].T,'P',second['cold_out'].p,second['cold_out'].fluid)
            second['hot_in'].q = (second['hot_out'].h - second['hot_in'].h)*second['hot_in'].m
            second['hot_out'].q = second['hot_in'].q
        elif process_case == 4:
            second['hot_in'].h = CP.PropsSI('H','T',second['hot_in'].T,'P',second['hot_in'].p,second['hot_in'].fluid)
            second['hot_out'].h = CP.PropsSI('H','T',second['hot_out'].T,'P',second['hot_out'].p,second['hot_out'].fluid)
            second['cold_in'].h = CP.PropsSI('H','T',second['cold_in'].T,'P',second['cold_in'].p,second['cold_in'].fluid)
            second['hot_in'].q = (second['hot_out'].h - second['hot_in'].h)*second['hot_in'].m
            second['hot_out'].q = second['hot_in'].q
        elif process_case == 5:
            second['hot_in'].h = CP.PropsSI('H','T',second['hot_in'].T,'P',second['hot_in'].p,second['hot_in'].fluid)
            second['hot_out'].h = CP.PropsSI('H','T',second['hot_out'].T,'P',second['hot_out'].p,second['hot_out'].fluid)
            second['cold_in'].h = CP.PropsSI('H','T',second['cold_in'].T,'P',second['cold_in'].p,second['cold_in'].fluid)
            second['cold_out'].h = CP.PropsSI('H','T',second['cold_out'].T,'P',second['cold_out'].p,second['cold_out'].fluid)
            second['hot_in'].q = (second['hot_out'].h - second['hot_in'].h)*second['hot_in'].m
            second['hot_out'].q = second['hot_in'].q
        
        return (second, process_case)
    
    def Solver(self, first, second, settings, outputs, max_p, process_case):
        PR_lb = 1.0
        PR_ub = 50.0
        beta = 0.00001
        a = 1
        results_array = []
        while a: 
            PR = 0.5*(PR_lb + PR_ub)
            for i in range(2):
                if i == 1:
                    PR = 0.5*(PR_lb + PR_ub)*(1+beta)
                    dPR = 0.5*(PR_lb + PR_ub)*beta

                first['hot_in'].p = max_p
                first['hot_out'].p = first['hot_in'].p*(1-settings.hot_dp)
                first['comp_in'].p = first['hot_in'].p / PR
                first['cold_out'].p = first['comp_in'].p /(1-settings.ihx_cold_dp)
                first['cold_in'].p = first['cold_out'].p / (1-settings.cold_dp)
                first['expand_in'].p = first['hot_out'].p*(1-settings.ihx_hot_dp)
                
                b = 1
                T_hout_lb = second['hot_in'].T
                T_hout_ub = min(second['hot_in'].T + 300.0, 1000)
                
                T_cout_lb = max(second['cold_in'].T - 300.0, 200.0)
                T_cout_ub = second['cold_in'].T 
                
                while b:
                    c = 1
                    while c:
                        T_hout = 0.5*(T_hout_lb+T_hout_ub)
                        first['hot_out'].T = T_hout
                        first['hot_out'].h = CP.PropsSI('H','T',first['hot_out'].T,'P',first['hot_out'].p,first['hot_out'].fluid)
                        T_cout = 0.5*(T_cout_lb+T_cout_ub)
                        first['cold_out'].T = T_cout
                        first['cold_out'].h = CP.PropsSI('H','T',first['cold_out'].T,'P',first['cold_out'].p,first['cold_out'].fluid)
                    
                        ideal_dh = min(first['hot_out'].h - CP.PropsSI('H','T',first['cold_out'].T,'P',first['expand_in'].p,first['hot_out'].fluid),CP.PropsSI('H','T',first['hot_out'].T,'P',first['comp_in'].p,first['cold_out'].fluid)-first['cold_out'].h)
                        first['comp_in'].h = first['cold_out'].h + settings.ihx_eff*ideal_dh
                        first['comp_in'].T = CP.PropsSI('T','H',first['comp_in'].h,'P',first['comp_in'].p,first['comp_in'].fluid)
                        first['comp_in'].s = CP.PropsSI('S','H',first['comp_in'].h,'P',first['comp_in'].p,first['comp_in'].fluid)
                        comp_out_h_ideal = CP.PropsSI('H','P',first['hot_in'].p,'S',first['comp_in'].s,first['hot_in'].fluid)
                        outputs.comp_W = (comp_out_h_ideal - first['comp_in'].h)/settings.comp_eff
                        
                        first['hot_in'].h = first['comp_in'].h + outputs.comp_W
                        first['hot_in'].T = CP.PropsSI('T','H',first['hot_in'].h, 'P',first['hot_in'].p, first['hot_in'].fluid)
                        
                        first['expand_in'].h = first['hot_out'].h - settings.ihx_eff*ideal_dh
                        first['expand_in'].T = CP.PropsSI('T','H',first['expand_in'].h,'P',first['expand_in'].p,first['expand_in'].fluid)
                        first['expand_in'].s = CP.PropsSI('S','H',first['expand_in'].h,'P',first['expand_in'].p,first['expand_in'].fluid)
                        expand_out_h_ideal = CP.PropsSI('H','P',first['cold_in'].p,'S',first['expand_in'].s,first['cold_in'].fluid)
                        outputs.expand_W = (first['expand_in'].h - expand_out_h_ideal)*settings.expand_eff
                        
                        first['cold_in'].h = first['expand_in'].h - outputs.expand_W
                        first['cold_in'].T = CP.PropsSI('T','H',first['cold_in'].h, 'P',first['cold_in'].p, first['cold_in'].fluid)
                        
                        if process_case <= 2:
                            first['cold_in'].q = -second['cold_in'].q
                            first['cold_in'].m = first['cold_in'].q / (first['cold_out'].h - first['cold_in'].h)
                            first['cold_out'].m = first['cold_in'].m
                            first['hot_in'].m = first['cold_in'].m
                            first['hot_out'].m = first['cold_in'].m
                            (dT_cold, T_rvs_cold, outputs.T_cold_sec) = self.HX_module(first['cold_in'], first['cold_out'], second['cold_in'], second['cold_out'], settings)
                            
                            if T_rvs_cold == 1:
                                T_cout_ub = T_cout
                            else:    
                                T_cold_err = (settings.cold_T_pp - dT_cold)/settings.cold_T_pp
                                if T_cold_err < 0:
                                    T_cout_lb = T_cout
                                else:
                                    T_cout_ub = T_cout
                                
                                if abs(T_cold_err) < 1.0e-3:
                                    c = 0
                                else:
                                    if T_cout_ub - T_cout_lb <= 0.01:
                                        c = 0
                        else:
                            first['hot_in'].q = -second['hot_in'].q
                            first['hot_in'].m = first['hot_in'].q / (first['hot_out'].h - first['hot_in'].h)
                            first['hot_out'].m = first['hot_in'].m
                            first['cold_in'].m = first['hot_in'].m
                            first['cold_out'].m = first['hot_in'].m
                            (dT_hot, T_rvs_hot, outputs.T_hot_sec) = self.HX_module(first['hot_in'], first['hot_out'], second['hot_in'], second['hot_out'], settings)
                            
                            if T_rvs_hot == 1:
                                T_hout_lb = T_hout
                            else:    
                                T_hot_err = (settings.hot_T_pp - dT_hot)/settings.hot_T_pp
                                if T_hot_err < 0:
                                    T_hout_ub = T_hout
                                else:
                                    T_hout_lb = T_hout
                                
                                if abs(T_hot_err) < 1.0e-3:
                                    c = 0
                                else:
                                    if T_hout_ub - T_hout_lb <= 0.01:
                                        c = 0
                    
                    if process_case <= 2:
                        first['hot_in'].q = first['hot_in'].m*(first['hot_out'].h - first['hot_in'].h)
                        second['hot_in'].q = -first['hot_in'].q
                        if process_case == 0:
                            second['hot_in'].h = second['hot_out'].h - second['hot_in'].q/second['hot_out'].m
                            second['hot_in'].T = CP.PropsSI('T','H',second['hot_in'].h,'P',second['hot_in'].p,second['hot_in'].fluid)
                        elif process_case == 1:
                            second['hot_out'].h = second['hot_in'].h + second['hot_in'].q/second['hot_in'].m
                            second['hot_out'].T = CP.PropsSI('T','H',second['hot_out'].h,'P',second['hot_out'].p,second['hot_out'].fluid)
                        elif process_case == 2:
                            second['hot_in'].m = second['hot_in'].q/(second['hot_out'].h - second['hot_in'].h)
                        
                        (dT_hot, T_rvs_hot, outputs.T_hot_sec) = self.HX_module(first['hot_in'], first['hot_out'], second['hot_in'], second['hot_out'], settings)
                        
                        if T_rvs_hot == 1:
                            T_hout_lb = T_hout
                        else:    
                            T_hot_err = (settings.hot_T_pp - dT_hot)/settings.hot_T_pp
                            if T_hot_err < 0:
                                T_hout_ub = T_hout
                            else:
                                T_hout_lb = T_hout
                            
                            if abs(T_hot_err) < 1.0e-3:
                                b = 0
                            else:
                                if T_hout_ub - T_hout_lb <= 0.01:
                                    b = 0
                    else:
                        first['cold_in'].q = first['cold_in'].m*(first['cold_out'].h - first['cold_in'].h)
                        second['cold_in'].q = -first['cold_in'].q
                        if process_case == 3:
                            second['cold_in'].h = second['cold_out'].h + outputs.cold_Q/second['cold_in'].m
                            second['cold_in'].T = CP.PropsSI('T','H',second['cold_in'].h,'P',second['cold_in'].p,second['cold_in'].fluid)
                        elif process_case == 4:
                            second['cold_out'].h = second['cold_in'].h - outputs.cold_Q/second['cold_out'].m
                            second['cold_out'].T = CP.PropsSI('T','H',second['cold_out'].h,'P',second['cold_out'].p,second['cold_out'].fluid)
                        elif process_case == 5:
                            second['cold_in'].m = second['cold_in'].q/(second['cold_out'].h - second['cold_in'].h)
                        
                        first['cold_out'].m = first['cold_in'].m
                        (dT_cold, T_rvs_cold, outputs.T_cold_sec) = self.HX_module(first['cold_in'], first['cold_out'], second['cold_in'], second['cold_out'], settings)
                        
                        if T_rvs_cold == 1:
                            T_cout_ub = T_cout
                        else:    
                            T_cold_err = (settings.cold_T_pp - dT_cold)/settings.cold_T_pp
                            if T_cold_err < 0:
                                T_cout_lb = T_cout
                            else:
                                T_cout_ub = T_cout
                            
                            if abs(T_cold_err) < 1.0e-3:
                                b = 0
                            else:
                                if T_cout_ub - T_cout_lb <= 0.01:
                                    b = 0
                    
                COP_cal = (first['hot_in'].h - first['hot_out'].h)/(outputs.comp_W/settings.mech_eff - outputs.expand_W*settings.gen_eff)
                
                if i == 0:
                    COP_0 = COP_cal
                else:
                    COP_1 = COP_cal
                    COP_avg = 0.5*(COP_1 + COP_0)
                    dCOP = ((COP_1 - COP_0)/COP_avg)/(dPR/(0.25*(PR_lb + PR_ub)*(2+beta)))
                    
            if dCOP > 0.0:
                PR_lb = PR
            else:
                PR_ub = PR
            
            results_array.append([COP_avg, dCOP])
            
            if len(results_array) > 2:
                if results_array[-2][0] > results_array[-1][0] and results_array[-2][0] > results_array[-3][0]:
                    if (PR_ub - PR_lb)/(0.25*(PR_lb + PR_ub)*(2+beta)) < self.tol:
                        outputs.COP = results_array[-2][0]
                        a = 0
                elif abs(results_array[-1][1]) < self.tol:
                    if (PR_ub - PR_lb)/(0.25*(PR_lb + PR_ub)*(2+beta)) < self.tol:
                        outputs.COP = results_array[-1][0]
                        a = 0
            
        outputs.comp_W = outputs.comp_W*first['cold_in'].m
        outputs.expand_W = outputs.expand_W*first['hot_in'].m
        outputs.ihx_Q = (first['hot_out'].h - first['expand_in'].h)*first['hot_in'].m
        outputs.hot_Q = (first['hot_in'].h - first['hot_out'].h)*first['hot_in'].m
        outputs.cold_Q = (first['cold_out'].h - first['cold_in'].h)*first['cold_in'].m
        
        hot_h_out_ideal = CP.PropsSI('H','T',second['hot_in'].T,'P',first['hot_out'].p,first['hot_out'].fluid)
        outputs.hot_eff = outputs.hot_Q/(first['hot_in'].h - hot_h_out_ideal)/first['hot_in'].m
        
        cold_h_out_ideal = CP.PropsSI('H','T',second['cold_in'].T,'P',first['cold_out'].p,first['cold_out'].fluid)
        outputs.cold_eff = outputs.cold_Q/(cold_h_out_ideal - first['cold_in'].h )/first['cold_in'].m
        
        return (first, second, outputs)
    
    def HX_module(self, first_in, first_out, second_in, second_out, inputs):
        h_primary = np.zeros(shape=(inputs.N_element+1))
        T_primary = np.zeros(shape=(inputs.N_element+1))
        P_primary = np.zeros(shape=(inputs.N_element+1))
        T_secondary = np.zeros(shape=(inputs.N_element+1))
        h_secondary = np.zeros(shape=(inputs.N_element+1))
        P_secondary = np.zeros(shape=(inputs.N_element+1))
        dT = np.zeros(shape=(inputs.N_element+1))
        
        P_primary[0]= first_in.p
        T_primary[0]= first_in.T
        h_primary[0]= first_in.h
        h_secondary[0]= second_out.h
        T_secondary[0]= second_out.T
        P_secondary[0]= second_out.p
        
        dT[0] = T_primary[0] - T_secondary[0]
        
        if (dT[0] < 0.0 and first_in.q < 0.0) or (dT[0] > 0.0 and first_in.q > 0.0):
            T_rvs = 1
        else:
            T_rvs = 0
            
        for i in range(inputs.N_element):
            if T_rvs == 1:
                break
            
            h_primary[i+1] = h_primary[i] + first_in.q/inputs.N_element/first_in.m
            P_primary[i+1] = P_primary[i] - (first_in.p - first_out.p)/inputs.N_element
            T_primary[i+1] = CP.PropsSI("T","H",h_primary[i+1],"P",P_primary[i+1],first_in.fluid)
            
            h_secondary[i+1] = h_secondary[i] + first_in.q/inputs.N_element/second_in.m
            P_secondary[i+1] = P_secondary[i] + (second_in.p - second_out.p)/inputs.N_element
            T_secondary[i+1] = CP.PropsSI("T","H",h_secondary[i+1],"P",P_secondary[i+1],second_in.fluid)
            
            dT[i+1] = T_primary[i+1] - T_secondary[i+1]
            
            if (dT[i+1] < 0.0 and first_in.q < 0.0) or (dT[i+1] > 0.0 and first_in.q > 0.0):
                T_rvs = 1
            else:
                T_rvs = 0
                
        T_pp = min(abs(dT))
        
        return(T_pp, T_rvs, T_secondary)
    
    def Plot_diagram(self, first, outputs):
        
        key_list = ['cold_out', 'comp_in', 'hot_in', 'hot_out', 'expand_in', 'cold_in', 'cold_out']
        
        T_vec = [first[t].T-273.15 for t in key_list]
        P_vec = [first[p].p/1.0e3 for p in key_list]
        for s in key_list:
            first[s].s = CP.PropsSI("S","T",first[s].T, "P", first[s].p, first['cold_out'].fluid)
        
        s_vec = [first[s].s/1.0e3 for s in key_list]
        h_vec = [first[h].h/1.0e3 for h in key_list]
        
        T_cold_sec = [t-273.15 for t in outputs.T_cold_sec]
        T_hot_sec = [t-273.15 for t in outputs.T_hot_sec]
        s_cold_sec = [(first['cold_in'].s+(first['cold_out'].s - first['cold_in'].s)*s/len(outputs.T_cold_sec))/1.0e3 for s in range(len(outputs.T_cold_sec))]
        s_hot_sec = [(first['hot_out'].s+(first['hot_in'].s - first['hot_out'].s)*s/len(outputs.T_hot_sec))/1.0e3 for s in range(len(outputs.T_hot_sec))]
        
        fig_ph, ax_ph = plt.subplots()
        ax_ph.plot(h_vec, P_vec, 'ko-')
        ax_ph.set_xlabel('Enthalpy [kJ/kg]',fontsize = 15)
        ax_ph.set_ylabel('Pressure [kPa]',fontsize = 15)
        ax_ph.set_title('Pressure-Enthalpy Diagram (Air)')
        ax_ph.tick_params(axis = 'x', labelsize = 13)
        ax_ph.tick_params(axis = 'y', labelsize = 13)
        
        fig_ts, ax_ts = plt.subplots()
        ax_ts.plot(s_vec, T_vec, 'ko-')
        ax_ts.plot(s_cold_sec, T_cold_sec, 'b', s_hot_sec, T_hot_sec[::-1], 'r')
        ax_ts.set_xlabel('Entropy [kJ/kg-K]',fontsize = 15)
        ax_ts.set_ylabel('Temperature [℃]',fontsize = 15)
        ax_ts.set_title('Temperature-Entropy Diagram (Air)')
        ax_ts.tick_params(axis = 'x', labelsize = 13)
        ax_ts.tick_params(axis = 'y', labelsize = 13)
        
        fig_ph.savefig('.\Ph_diagram.png',dpi=300)
        fig_ts.savefig('.\Ts_diagram.png',dpi=300)
        
    def Post_process(self, first, second, settings, outputs):
        print("COP: %.3f" %(outputs.COP))
        print("Q_hot: %.3f [kW]" %(outputs.hot_Q/1.0e3))
        print("Q_cold: %.3f [kW]" %(outputs.cold_Q/1.0e3))
        print("W_comp: %.3f [kW]" %(outputs.comp_W/settings.mech_eff/1.0e3))
        print("W_expand: %.3f [kW]" %(outputs.expand_W*settings.gen_eff/1.0e3))
        print("Comp_eff: %.3f [%%], Expand_eff: %.3f [%%], Mech_eff: %.3f [%%], Gen_eff: %.3f [%%]" %(settings.comp_eff*100, settings.expand_eff*100, settings.mech_eff*100, settings.gen_eff*100))
        print("Compression ratio: %.3f" %(first['hot_in'].p/first['cold_out'].p))
        print("[First side cold state] T_in: %.3f [℃], P_in: %.3f [kPa], T_out: %.3f [℃], P_out: %.3f [kPa]" %(first['cold_in'].T-273.15, first['cold_in'].p/1.0e3, first['cold_out'].T-273.15, first['cold_out'].p/1.0e3))
        print("[First side hot state] T_in: %.3f [℃], P_in: %.3f [kPa], T_out: %.3f [℃], P_out: %.3f [kPa], mdot: %.3f [kg/s]" %(first['hot_in'].T-273.15, first['hot_in'].p/1.0e3, first['hot_out'].T-273.15, first['hot_out'].p/1.0e3, first['cold_in'].m))
        print("[Process side cold state] T_in: %.3f [℃], P_in: %.3f [kPa], T_out: %.3f [℃], P_out: %.3f [kPa], mdot: %.3f [kg/s]" %(second['cold_in'].T-273.15, second['cold_in'].p/1.0e3, second['cold_out'].T-273.15, second['cold_out'].p/1.0e3, second['cold_in'].m))
        print("[Process side hot state] T_in: %.3f [℃], P_in: %.3f [kPa], T_out: %.3f [℃], P_out: %.3f [kPa], mdot: %.3f [kg/s]" %(second['hot_in'].T-273.15, second['hot_in'].p/1.0e3, second['hot_out'].T-273.15, second['hot_out'].p/1.0e3, second['hot_in'].m))
        print('-------------------------------------------------------------------')
    
if __name__ == "__main__":
    
    max_p = 5.0e5
    
    hot_fluid = 'water'
    cold_fluid = 'Water'
    ref_fluid = 'air'
    
    hot_in_T = 350+273.15
    hot_in_p = 18.0e6
    hot_out_T = 360+273.15
    hot_out_p = 18.0e6
    hot_m = 1.0
    
    cold_in_T = 270+273.15
    cold_out_T = 260+273.15
    cold_in_p = 5.0e6
    cold_out_p = 5.0e6
    cold_m = 0.0
    
    hot_in_sec = Fluid_stream(T=hot_in_T, p=hot_in_p, m=hot_m, fluid=hot_fluid)
    hot_out_sec = Fluid_stream(T=hot_out_T, p=hot_out_p, m=hot_m, fluid=hot_fluid)
    cold_in_sec = Fluid_stream(T=cold_in_T, p=cold_in_p, m=cold_m, fluid=cold_fluid)
    cold_out_sec = Fluid_stream(T=cold_out_T, p=cold_out_p, m=cold_m, fluid=cold_fluid)
    hot_in_pri = Fluid_stream(fluid = ref_fluid) 
    hot_out_pri = Fluid_stream(fluid = ref_fluid)
    cold_in_pri = Fluid_stream(fluid = ref_fluid)
    cold_out_pri = Fluid_stream(fluid = ref_fluid)
    comp_in_pri = Fluid_stream(fluid = ref_fluid)
    expand_in_pri = Fluid_stream(fluid = ref_fluid)
    
    first = {'hot_in':hot_in_pri, 'hot_out':hot_out_pri, 'cold_in': cold_in_pri, 'cold_out':cold_out_pri, 'comp_in':comp_in_pri, 'expand_in':expand_in_pri}
    second = {'hot_in':hot_in_sec, 'hot_out':hot_out_sec, 'cold_in':cold_in_sec, 'cold_out':cold_out_sec}
    
    settings = Settings(hot_dp = 0.01, hot_T_pp=5, cold_dp = 0.01, cold_T_pp = 5.0, ihx_cold_dp = 0.01, ihx_hot_dp = 0.01, ihx_eff = 0.9, comp_eff = 0.75, expand_eff = 0.85, gen_eff=0.99, N_element=30)
    outputs = Outputs()
    
    airHP = Brayton()
    [second, process_case] = airHP.Input_Processing(second)
    airHP.Solver(first, second, settings, outputs, max_p, process_case)
    airHP.Plot_diagram(first, outputs)
    airHP.Post_process(first, second, settings, outputs)