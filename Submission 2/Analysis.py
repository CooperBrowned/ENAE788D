import math
import numpy as np

#givens
G0 = 9.81              
DV_TOTAL = 12300        
M_PAYLOAD = 26000         
PAYLOAD_D = 5.2              
PAYLOAD_H = 13   
LD_MAX = 13               
TW_MIN = {1: 1.3, 2: 0.76}   
MARGIN = 0.30 
ENGINE_LENGTH = 3  
FUEL_MIXTURE_MASS_RATIO = [3.6, 6.03, 2.72, 0.0, 2.67]
THRUST_PER_MOTOR_STAGE_1 = [2.26, 1.86, 1.92, 4.5, 1.75]
THRUST_PER_MOTOR_STAGE_2 = [0.745, 0.099, 0.061, 2.94, 0.067]
EXHAUST_DIAMETER_STAGE_1 = [2.4, 2.4, 3.7, 6.6, 1.5]
EXHAUST_DIAMETER_STAGE_2 = [1.5, 2.15, 0.92, 2.34, 1.13]
CHAMBER_PRESSURE_STAGE_1 = [35.16, 20.64, 25.8, 10.5, 15.7]
CHAMBER_PRESSURE_STAGE_2 = [10.1, 4.2, 6.77, 5.0, 14.7]
NOZZLE_EXPANSION_STAGE_1 = [34.34, 78.0, 37.0, 16.0, 26.2]
NOZZLE_EXPANSION_STAGE_2 = [45.0, 84.0, 14.5, 56.0, 81.3]
SIGMA = [0.067, 0.075, 0.063, 0.087, 0.061]
ISP_SEA = [327.0, 366.0, 311.0, 269.0, 285.0]
ISP_VAC = [380.0, 452.0, 337.0, 279.0, 316.0]
DENSITY = {"LH2": 71.0, "LOX": 1140.0, "RP1": 820.0, "LCH4": 423.0,"APCP": 1680.0, "N2O4": 1442.0, "UDMH": 791.0}
PROPELLANT_PAIRS = [("LCH4", "LOX"), ("LH2", "LOX"), ("RP1", "LOX"), ("APCP", None), ("UDMH", "N2O4")]

#assumptions
DIAMETER = [5.2,5.2]
INTERSTAGE_LENGTH = 3  
AFT_FAIRING_LENGTH = 3 
# assuming wiring length is stage 1: tank + aft fairing; stage 2: tank + engine section + payload fairing (no double counting)
# assuming both aft and interstage fairings stay attached to first stage
# assuming payload fairing is not jettisoned until after 2nd stage burn finishes
# assuming tanks are spherical unless diameter is greater than rocket diameter, at which point they become perfect hemispherical

NOSECONE_HEIGHT = PAYLOAD_D
PAYLOAD_FAIRING_LENGTH = PAYLOAD_H + NOSECONE_HEIGHT  # cylinder + nosecone (m)

payload_fairing_mass = 4.95*((PAYLOAD_D*math.pi*PAYLOAD_H) + (math.pi*PAYLOAD_D/2)*math.sqrt((PAYLOAD_D/2)**2+(NOSECONE_HEIGHT)**2))**1.15

def update_stage(current_stage_mass, deltaV, propellant_index, stage, payload_mass):
    propellant_mass = get_propellant_mass(stage,propellant_index,current_stage_mass,deltaV)
    thrust_required = current_stage_mass*G0*TW_MIN[stage]
    tank_mass, tank_length = get_tank_mass(propellant_index, propellant_mass, stage)
    insulation_mass = get_insulation_mass(propellant_index, propellant_mass, stage)
    casing_mass = get_casing_mass(propellant_index, propellant_mass, stage)
    engine_mass, n_engines = get_engine_mass(propellant_index, thrust_required, stage)
    gimbal_mass = get_gimbal_mass(propellant_index, thrust_required, stage)
    thrust_structure_mass =2.55*10**-4*thrust_required
    avionics_mass = 10*current_stage_mass**0.361
    stage_payload_fairing_mass = 0.0 if stage == 1 else payload_fairing_mass  # fairing is carried by stage 2
    stage_length = get_stage_length(tank_length, stage)
    wiring_mass = 1.058*math.sqrt(current_stage_mass)*stage_length**0.25
    intertank_mass = get_intertank_fairing_mass(propellant_index, propellant_mass, stage)
    interstage_mass = get_interstage_mass(stage)
    aft_fairing_mass = get_aft_fairing_mass(stage)
    inert_mass = (1 + MARGIN)*(tank_mass + insulation_mass + casing_mass + engine_mass + gimbal_mass + thrust_structure_mass + avionics_mass + wiring_mass + stage_payload_fairing_mass+ intertank_mass + interstage_mass + aft_fairing_mass)
    total_stage_mass = inert_mass + propellant_mass
    stage_vars = {
        "tank_mass": tank_mass,
        "insulation_mass": insulation_mass,
        "casing_mass": casing_mass,
        "engine_mass": engine_mass,
        "gimbal_mass": gimbal_mass,
        "thrust_structure_mass": thrust_structure_mass,
        "avionics_mass": avionics_mass,
        "wiring_mass": wiring_mass,
        "intertank_mass": intertank_mass,
        "interstage_mass": interstage_mass,
        "aft_fairing_mass": aft_fairing_mass,
        "n_engines": n_engines,
        "stage_length": stage_length,
        "propellant_mass": propellant_mass,
    }
    return total_stage_mass, stage_vars


def get_stage_length(tank_length, stage):
    # length used for the wiring MER; each section is counted once, so upper stages are not included
    if stage == 1:
        return tank_length + AFT_FAIRING_LENGTH
    return tank_length + ENGINE_LENGTH + PAYLOAD_FAIRING_LENGTH

def get_component_masses(propellant_index, propellant_mass):
    fuel, ox = PROPELLANT_PAIRS[propellant_index]
    of_ratio = FUEL_MIXTURE_MASS_RATIO[propellant_index]
    if ox is None:
        return {fuel: propellant_mass}
    return {fuel: propellant_mass/(1 + of_ratio), ox: propellant_mass*of_ratio/(1 + of_ratio)}

def get_stage_diameter(stage):
    return DIAMETER[stage - 1]

def get_fairing_mass(stage, length):
    return 4.95*(math.pi*get_stage_diameter(stage)*length)**1.15

def get_interstage_mass(stage):
    # both fairings are carried by stage 1 only, so stage 2 has none
    if stage != 1:
        return 0.0
    return get_fairing_mass(stage, INTERSTAGE_LENGTH)

def get_aft_fairing_mass(stage):
    if stage != 1:
        return 0.0
    return get_fairing_mass(stage, AFT_FAIRING_LENGTH)

def get_tank_geometry(volume, stage):
    max_radius = get_stage_diameter(stage)/2
    sphere_radius = (3*volume/(4*math.pi))**(1/3)
    if sphere_radius <= max_radius:
        return 2*sphere_radius, 4*math.pi*sphere_radius**2, sphere_radius
    cylinder_length = (volume - 4/3*math.pi*max_radius**3)/(math.pi*max_radius**2)
    length = cylinder_length + 2*max_radius
    area = 2*math.pi*max_radius*cylinder_length + 4*math.pi*max_radius**2
    return length, area, max_radius

def get_intertank_fairing_mass(propellant_index, propellant_mass, stage):
    if is_solid(propellant_index):
        return 0.0
    cap_heights = sum(get_tank_geometry(mass/DENSITY[name], stage)[2]
                      for name, mass in get_component_masses(propellant_index, propellant_mass).items())
    return get_fairing_mass(stage, cap_heights)

def is_solid(propellant_index):
    #Simple return to clean up if statements for weird solid edge cases
    return PROPELLANT_PAIRS[propellant_index][1] is None

def get_casing_mass(propellant_index, propellant_mass, stage):
    if is_solid(propellant_index):
        return 0.135*propellant_mass
    else:
        return 0.0

def get_tank_mass(propellant_index, propellant_mass, stage):
    tank_mass = 0.0
    tank_length = 0.0
    if is_solid(propellant_index):
        volume = propellant_mass/DENSITY[PROPELLANT_PAIRS[propellant_index][0]]
        tank_length += get_tank_geometry(volume, stage)[0]
        return tank_mass, tank_length
    for name, mass in get_component_masses(propellant_index, propellant_mass).items():
        volume = mass/DENSITY[name]
        tank_mass += (9.09 if name == "LH2" else 12.16)*volume
        tank_length += get_tank_geometry(volume, stage)[0]
    return tank_mass, tank_length

def get_insulation_mass(propellant_index, propellant_mass, stage):
    insulation_density = {"LH2": 2.88, "LOX": 1.123, "LCH4": 1.123} 
    insulation_mass = 0.0
    for name, mass in get_component_masses(propellant_index, propellant_mass).items():
        if name in insulation_density:
            insulation_mass += insulation_density[name]*get_tank_geometry(mass/DENSITY[name], stage)[1]
    return insulation_mass
    
def get_engine_mass(propellant_index, thrust_required, stage):
    if stage == 1:
        thrust_table = THRUST_PER_MOTOR_STAGE_1 
    else:
        thrust_table = THRUST_PER_MOTOR_STAGE_2
    if stage == 1:
        expansion_table = NOZZLE_EXPANSION_STAGE_1  
    else:
        expansion_table =  NOZZLE_EXPANSION_STAGE_2
    thrust_per_motor = thrust_table[propellant_index]*1e6  
    expansion_ratio = expansion_table[propellant_index]
    n_engines = math.ceil(thrust_required/thrust_per_motor)
    mass_per_engine = (7.81e-4*thrust_per_motor+ 3.37e-5*thrust_per_motor*math.sqrt(expansion_ratio)+ 59)
    if is_solid(propellant_index):
        return 0.0, n_engines
    return n_engines*mass_per_engine, n_engines

def get_gimbal_mass(propellant_index, thrust_required, stage):
    if stage == 1:
        thrust_per_motor = THRUST_PER_MOTOR_STAGE_1[propellant_index]*1e6  # MN -> N
        chamber_pressure = CHAMBER_PRESSURE_STAGE_1[propellant_index]*1e6  # MPa -> Pa
    else:
        thrust_per_motor = THRUST_PER_MOTOR_STAGE_2[propellant_index]*1e6
        chamber_pressure = CHAMBER_PRESSURE_STAGE_2[propellant_index]*1e6
    n_motors = math.ceil(thrust_required/thrust_per_motor)
    return n_motors*237.8*(thrust_per_motor/chamber_pressure)**0.9375

def get_propellant_mass(stage, propellant_index, current_stage_mass, deltaV):
    if stage == 1:
        isp = ISP_SEA[propellant_index]
    else:
        isp = ISP_VAC[propellant_index]
    return(current_stage_mass*(1 - math.exp(-deltaV / (isp*G0))))

#adapted code from matlab explicit solution
# def initialize(delta_v1, stage1_idx, stage2_idx):
#     isp_sea = ISP_SEA[stage1_idx]
#     isp_vac = ISP_VAC[stage2_idx]
#     sigma1 = SIGMA[stage1_idx]
#     sigma2 = SIGMA[stage2_idx]
#     e1 = math.exp(delta_v1 / (G0 * isp_sea))
#     e2 = math.exp((delta_v1 - DV_TOTAL) / (G0 * isp_vac))
#     denom = (sigma1 * e1 - 1.0) * (sigma2 - e2)
#     m1 = (M_PAYLOAD * (e1 - 1.0) + M_PAYLOAD * sigma1 * e1) / denom
#     m2 = -(M_PAYLOAD - M_PAYLOAD * e2) / (sigma2 - e2) - (M_PAYLOAD * sigma2) / (sigma2 - e2)
#     return [m1, m2]


def converge_stage(above_mass, deltaV, propellant_index, stage, tol=1e-6, max_iter=200): 
    m0 = 1  
    for _ in range(max_iter):
        new_stage_mass, stage_vars = update_stage(m0, deltaV, propellant_index, stage, above_mass)
        new_m0 = new_stage_mass + above_mass
        if abs(new_m0 - m0) < tol*m0:
            return new_stage_mass, stage_vars
        m0 = new_m0
    raise RuntimeError("stage did not converge") 


def size_vehicle(delta_v1, stage1_idx, stage2_idx):
    delta_v2 = DV_TOTAL - delta_v1
    #guess1, guess2 = initialize(delta_v1, stage1_idx, stage2_idx)   
    stage2_mass, stage2_vars = converge_stage(M_PAYLOAD, delta_v2, stage2_idx, 2) 
    stage1_mass, stage1_vars = converge_stage(stage2_mass + M_PAYLOAD, delta_v1, stage1_idx, 1)
    if (stage1_vars["stage_length"]+stage2_vars["stage_length"])/DIAMETER[0] > LD_MAX:
        raise RuntimeError("stage L/D exceeds limit"+" stage1_length: "+str(stage1_vars["stage_length"])+" stage2_length: "+str(stage2_vars["stage_length"]))
    if (stage1_vars["n_engines"])*(EXHAUST_DIAMETER_STAGE_1[stage1_idx]/2)**2*math.pi > (DIAMETER[0]/2)**2*math.pi:
        raise RuntimeError("stage 1 engine diameter exceeds limit")
    if (stage2_vars["n_engines"])*(EXHAUST_DIAMETER_STAGE_2[stage2_idx]/2)**2*math.pi > (DIAMETER[1]/2)**2*math.pi:
        raise RuntimeError("stage 2 engine diameter exceeds limit")
    cost = 13.52 * (stage1_mass - stage1_vars["propellant_mass"])**0.55 + 13.52 * (stage2_mass - stage2_vars["propellant_mass"])**0.55
    return stage1_mass, stage2_mass, stage1_mass + stage2_mass + M_PAYLOAD, cost, stage1_vars, stage2_vars


def print_vehicle(diameter, result):
    stage1_mass, stage2_mass, total_mass, cost, stage1_vars, stage2_vars = result
    print("diameter:", diameter)
    for name in stage1_vars:
        print(name, round(stage1_vars[name], 1), round(stage2_vars[name], 1))
    print("stage masses:", round(stage1_mass), round(stage2_mass))
    print("total mass:", round(total_mass), "cost:", round(cost))


if __name__ == "__main__":
    stage_1_idx = 1
    stage_2_idx = 0
    #deltaV1 = 5720
    deltaV1 = 6396
    results = []
    #itterate from deltaV = 1000 to 10000
    
    for DIAMETER[0] in np.linspace(5.2, 15, 100):
        DIAMETER[1] = DIAMETER[0]
        try:
            results.append([DIAMETER[0],size_vehicle(deltaV1, stage_1_idx, stage_2_idx)])
        except RuntimeError as e:
            print(e)
    #finds the minimum overall mass with the corresponding diameter
    min_mass_d, min_mass_r = min(results, key=lambda r: r[1][2])
    print("\n\n\nMinimum mass results:")
    print_vehicle(min_mass_d, min_mass_r)
    min_cost_d, min_cost_r = min(results, key=lambda r: r[1][3])
    print("\n\n\nMinimum cost results:")
    print_vehicle(min_cost_d, min_cost_r)