import math

G0 = 9.81              
DV_TOTAL = 12300        
M_PAYLOAD = 26000         
PAYLOAD_D = 5.2              
PAYLOAD_H = 13.0             
LD_MAX = 13.0                
TW_MIN = {1: 1.3, 2: 0.76}   

MARGIN = 0.30      

FUEL_MIXTURE_MASS_RATIO = [3.6, 6.03, 2.72, 0.0, 2.67]
THRUST_PER_MOTOR_STAGE_1 = [2.26, 1.86, 1.92, 4.5, 1.75]
THRUST_PER_MOTOR_STAGE_2 = [0.745, 0.099, 0.061, 2.94, 0.067]
EXHAUST_DIAMETER_STAGE_1 = [2.4, 2.4, 3.7, 6.6, 1.5]
EXHAUST_DIAMETER_STAGE_2 = [1.5, 2.15, 0.92, 2.34, 1.13]
CHAMBER_PRESSURE_STAGE_1 = [35.16, 20.64, 25.8, 10.5, 15.7]
CHAMBER_PRESSURE_STAGE_2 = [10.1, 4.2, 6.77, 5.0, 14.7]
NOZZLE_EXPANSION_STAGE_1 = [34.34, 78.0, 37.0, 16.0, 26.2]
NOZZLE_EXPANSION_STAGE_2 = [45.0, 84.0, 14.5, 56.0, 81.3]
RHO_PROPELLANT = [71.0, 1140.0, 820.0, 423.0, 1680.0, 1442.0, 791.0]
SIGMA = [0.067, 0.075, 0.063, 0.087, 0.061]
ISP_SEA = [327.0, 366.0, 311.0, 269.0, 285.0]
ISP_VAC = [380.0, 452.0, 337.0, 279.0, 316.0]
          
 
DENSITY = {"LH2": 71.0, "LOX": 1140.0, "RP1": 820.0, "LCH4": 423.0,"APCP": 1680.0, "N2O4": 1442.0, "UDMH": 791.0}





#adapted code from matlab explicit solution
def initialize(delta_v1, stage1_idx, stage2_idx):
    isp_sea = ISP_SEA[stage1_idx]
    isp_vac = ISP_VAC[stage2_idx]
    sigma1 = SIGMA[stage1_idx]
    sigma2 = SIGMA[stage2_idx]
    e1 = math.exp(delta_v1 / (G0 * isp_sea))
    e2 = math.exp((delta_v1 - DV_TOTAL) / (G0 * isp_vac))
    denom = (sigma1 * e1 - 1.0) * (sigma2 - e2)
    m1 = (M_PAYLOAD * (e1 - 1.0) + M_PAYLOAD * sigma1 * e1) / denom
    m2 = -(M_PAYLOAD - M_PAYLOAD * e2) / (sigma2 - e2) - (M_PAYLOAD * sigma2) / (sigma2 - e2)
    return [m1, m2]
