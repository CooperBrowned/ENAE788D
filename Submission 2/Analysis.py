import math
import numpy as np

# givens
G0 = 9.81                 # m/s^2
DV_TOTAL = 12300          # m/s
M_PAYLOAD = 26000         # kg
PAYLOAD_D = 5.2           # m
PAYLOAD_H = 13            # m
LD_MAX = 13
TW_MIN = {1: 1.3, 2: 0.76}
MARGIN = 0.30
ENGINE_LENGTH = 3         # m

# propellant tables, indexed by propellant_index (see PROPELLANT_PAIRS)
FUEL_MIXTURE_MASS_RATIO = [3.6, 6.03, 2.72, 0.0, 2.67]   # O/F
THRUST_PER_MOTOR_STAGE_1 = [2.26, 1.86, 1.92, 4.5, 1.75]  # MN
THRUST_PER_MOTOR_STAGE_2 = [0.745, 0.099, 0.061, 2.94, 0.067]
EXHAUST_DIAMETER_STAGE_1 = [2.4, 2.4, 3.7, 6.6, 1.5]  # m
EXHAUST_DIAMETER_STAGE_2 = [1.5, 2.15, 0.92, 2.34, 1.13]
EXHAUST_DIAMETER = [EXHAUST_DIAMETER_STAGE_1, EXHAUST_DIAMETER_STAGE_2]
CHAMBER_PRESSURE_STAGE_1 = [35.16, 20.64, 25.8, 10.5, 15.7]  # MPa
CHAMBER_PRESSURE_STAGE_2 = [10.1, 4.2, 6.77, 5.0, 14.7]
NOZZLE_EXPANSION_STAGE_1 = [34.34, 78.0, 37.0, 16.0, 26.2]
NOZZLE_EXPANSION_STAGE_2 = [45.0, 84.0, 14.5, 56.0, 81.3]
SIGMA = [0.067, 0.075, 0.063, 0.087, 0.061]
ISP_SEA = [327.0, 366.0, 311.0, 269.0, 285.0]  # s
ISP_VAC = [380.0, 452.0, 337.0, 279.0, 316.0]
DENSITY = {"LH2": 71.0, "LOX": 1140.0, "RP1": 820.0, "LCH4": 423.0, "APCP": 1680.0, "N2O4": 1442.0, "UDMH": 791.0}  # kg/m^3
PROPELLANT_PAIRS = [("LCH4", "LOX"), ("LH2", "LOX"), ("RP1", "LOX"), ("APCP", None), ("UDMH", "N2O4")]
ENGINE_PACKING_FACTOR = {1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
                         8: 3.304, 9: 3.613, 10: 3.813, 11: 3.923, 12: 4.03, 13: 4.236}
# assumptions
DIAMETER = [5.2, 5.2]
AFT_FAIRING_LENGTH = 3
# assuming wiring length is stage 1: tank + aft fairing; stage 2: tank + engine section + payload fairing (no double counting)
# assuming both aft and interstage fairings stay attached to first stage
# assuming payload fairing is not jettisoned until after 2nd stage burn finishes
# assuming tanks are spherical unless diameter is greater than rocket diameter, at which point they become a cylinder with hemispherical caps

NOSECONE_HEIGHT = PAYLOAD_D
PAYLOAD_FAIRING_LENGTH = PAYLOAD_H + NOSECONE_HEIGHT  # cylinder + nosecone (m)

# fairing MER applied to the cylinder area plus the cone's lateral area
payload_fairing_mass = 4.95*((PAYLOAD_D*math.pi*PAYLOAD_H) + (math.pi*PAYLOAD_D/2)*math.sqrt((PAYLOAD_D/2)**2 + NOSECONE_HEIGHT**2))**1.15


# FUNCTION CATEGORY: stage sizing
def update_stage(current_stage_mass, deltaV, propellant_index, stage):
    # Cooper, Ethan, Leela, Nihira, Colby
    # one pass of the stage sizing: given a guess for the stage's starting mass,
    # add up all the inert components and return the new stage mass
    propellant_mass = get_propellant_mass(stage, propellant_index, current_stage_mass, deltaV)
    thrust_required = current_stage_mass*G0*TW_MIN[stage]
    tank_mass, tank_length = get_tank_mass(propellant_index, propellant_mass, stage)
    insulation_mass = get_insulation_mass(propellant_index, propellant_mass, stage)
    casing_mass = get_casing_mass(propellant_index, propellant_mass, stage)
    engine_mass, n_engines = get_engine_mass(propellant_index, thrust_required, stage)
    gimbal_mass = get_gimbal_mass(propellant_index, thrust_required, stage)
    thrust_structure_mass = 2.55e-4*thrust_required
    avionics_mass = 10*current_stage_mass**0.361
    stage_payload_fairing_mass = 0.0 if stage == 1 else payload_fairing_mass  # fairing is carried by stage 2
    stage_length = get_stage_length(tank_length, stage)
    wiring_mass = 1.058*math.sqrt(current_stage_mass)*stage_length**0.25
    intertank_mass = get_intertank_fairing_mass(propellant_index, propellant_mass, stage)
    interstage_mass = get_interstage_mass(stage)
    aft_fairing_mass = get_aft_fairing_mass(stage)
    inert_mass = (1 + MARGIN)*(tank_mass + insulation_mass + casing_mass + engine_mass + gimbal_mass
                               + thrust_structure_mass + avionics_mass + wiring_mass + stage_payload_fairing_mass
                               + intertank_mass + interstage_mass + aft_fairing_mass)
    total_stage_mass = inert_mass + propellant_mass
    #pass stage vars out for later analysis / printing
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
        "payload_fairing_mass": stage_payload_fairing_mass,
    }
    return total_stage_mass, stage_vars

def converge_stage(above_mass, deltaV, propellant_index, stage, tol=1e-6, max_iter=2000):
    # Cooper
    #the stage mass depends on itself so iterate update_stage until the total liftoff mass stops changing or error out
    m0 = 1
    for _ in range(max_iter):
        new_stage_mass, stage_vars = update_stage(m0, deltaV, propellant_index, stage)
        new_m0 = new_stage_mass + above_mass
        if abs(new_m0 - m0) < tol*m0:
            return new_stage_mass, stage_vars
        m0 = new_m0
    raise RuntimeError("stage did not converge")


def size_vehicle(delta_v1, stage1_idx, stage2_idx):
    # Cooper
    #sizes stage 2 first, then stage 1 under it
    # throws if the vehicle braeks the L/D or engine fit constraints, otherwise returns masses and cost
    delta_v2 = DV_TOTAL - delta_v1
    stage2_mass, stage2_vars = converge_stage(M_PAYLOAD, delta_v2, stage2_idx, 2)
    stage1_mass, stage1_vars = converge_stage(stage2_mass + M_PAYLOAD, delta_v1, stage1_idx, 1)
    total_length = stage1_vars["stage_length"] + stage2_vars["stage_length"]
    if total_length/DIAMETER[0] > LD_MAX:
        raise RuntimeError(f"stage L/D exceeds limit stage1_length: {stage1_vars['stage_length']} "
                           f"stage2_length: {stage2_vars['stage_length']}")
    check_engines(stage1_vars["n_engines"], stage1_idx, 0)
    check_engines(stage2_vars["n_engines"], stage2_idx, 1)

    cost = 13.52*(stage1_mass - stage1_vars["propellant_mass"])**0.55 + 13.52*(stage2_mass - stage2_vars["propellant_mass"])**0.55
    return stage1_mass, stage2_mass, stage1_mass + stage2_mass + M_PAYLOAD, cost, stage1_vars, stage2_vars

# FUNCTION CATEGORY: stage property getters
def get_stage_length(tank_length, stage):
    # Cooper
    # length used for the wiring MER; each section is counted once, so upper stages are not included, from lecture
    if stage == 1:
        return tank_length + AFT_FAIRING_LENGTH
    return tank_length + ENGINE_LENGTH + PAYLOAD_FAIRING_LENGTH


def get_component_masses(propellant_index, propellant_mass):
    # Ethan
    # splits total propellant into fuel and oxidizer using the given ratio 
    fuel, ox = PROPELLANT_PAIRS[propellant_index]
    of_ratio = FUEL_MIXTURE_MASS_RATIO[propellant_index]
    if ox is None:
        return {fuel: propellant_mass}
    return {fuel: propellant_mass/(1 + of_ratio), ox: propellant_mass*of_ratio/(1 + of_ratio)}


def get_stage_diameter(stage):
    return DIAMETER[stage - 1]


def get_fairing_mass(stage, length):
    # fairing MER based on the cylinder at the stage diameter
    return 4.95*(math.pi*get_stage_diameter(stage)*length)**1.15


def get_interstage_mass(stage):
    # Nihira
    # both fairings are carried by stage 1 only, so stage 2 has none
    if stage != 1:
        return 0.0
    # interstage has to cover the stage 2 engine plus the bottom tank dome
    interstage_length = ENGINE_LENGTH + get_stage_diameter(stage)
    return get_fairing_mass(stage, interstage_length)


def get_aft_fairing_mass(stage):
    # Nihira
    # aft fairing covers the stage 1 engines, stage 2 is included in interstage
    if stage != 1:
        return 0.0
    return get_fairing_mass(stage, AFT_FAIRING_LENGTH)


def get_tank_geometry(volume, stage):
    # Ethan
    # returns (length, surface area, radius) of a tank holding the given volume.
    # sphere if it fits inside the stage, otherwise cylinder with hemispherical ends
    max_radius = get_stage_diameter(stage)/2
    sphere_radius = (3*volume/(4*math.pi))**(1/3)
    if sphere_radius <= max_radius:
        return 2*sphere_radius, 4*math.pi*sphere_radius**2, sphere_radius
    cylinder_length = (volume - 4/3*math.pi*max_radius**3)/(math.pi*max_radius**2)
    length = cylinder_length + 2*max_radius
    area = 2*math.pi*max_radius*cylinder_length + 4*math.pi*max_radius**2
    return length, area, max_radius


def get_intertank_fairing_mass(propellant_index, propellant_mass, stage):
    # Nihira
    # intertank spans the two facing tank domes, so length is the sum of the cap heights
    if is_solid(propellant_index):
        return 0.0
    cap_heights = sum(get_tank_geometry(mass/DENSITY[name], stage)[2]
                      for name, mass in get_component_masses(propellant_index, propellant_mass).items())
    return get_fairing_mass(stage, cap_heights)


def is_solid(propellant_index):
    # Leela
    # simple check to clean up if statements for the solid edge cases
    return PROPELLANT_PAIRS[propellant_index][1] is None


def get_casing_mass(propellant_index, propellant_mass, stage):
    # only solids have a motor casing, liquids use tanks instead
    if is_solid(propellant_index):
        return 0.135*propellant_mass
    return 0.0


def get_tank_mass(propellant_index, propellant_mass, stage):
    # Cooper
    # tank mass from volume (different MER for LH2) and the total tank stack length.
    # solids get no tank mass since the casing covers it, but we still need the length
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
    # Cooper
    # insulation only on the cryogenic tanks from lecture eqn
    insulation_density = {"LH2": 2.88, "LOX": 1.123, "LCH4": 1.123}
    insulation_mass = 0.0
    for name, mass in get_component_masses(propellant_index, propellant_mass).items():
        if name in insulation_density:
            insulation_mass += insulation_density[name]*get_tank_geometry(mass/DENSITY[name], stage)[1]
    return insulation_mass


def get_engine_mass(propellant_index, thrust_required, stage):
    # Cooper
    # picks the number of engines needed to hit the required thrust and returns (total mass, n_engines).
    # solid motor mass is already in the casing so it returns 0 mass but still counts motors
    #solids modeling in general is questionable here, but the givens are pretty clear in the problem statement. 
    #I would have liked to model them a single engine/nozzle given how they work, but oh well.
    if stage == 1:
        thrust_table = THRUST_PER_MOTOR_STAGE_1
        expansion_table = NOZZLE_EXPANSION_STAGE_1
    else:
        thrust_table = THRUST_PER_MOTOR_STAGE_2
        expansion_table = NOZZLE_EXPANSION_STAGE_2
    thrust_per_motor = thrust_table[propellant_index]*1e6
    expansion_ratio = expansion_table[propellant_index]
    #need the ceiling here because if you don't have enough thrust it doesn't work and no fracitonal engines
    n_engines = math.ceil(thrust_required/thrust_per_motor)
    mass_per_engine = 7.81e-4*thrust_per_motor + 3.37e-5*thrust_per_motor*math.sqrt(expansion_ratio) + 59
    if is_solid(propellant_index):
        return 0.0, n_engines
    return n_engines*mass_per_engine, n_engines


def get_gimbal_mass(propellant_index, thrust_required, stage):
    # Leela
    # gimbal MER per motor, depends on thrust and chamber pressure
    if stage == 1:
        thrust_per_motor = THRUST_PER_MOTOR_STAGE_1[propellant_index]*1e6
        chamber_pressure = CHAMBER_PRESSURE_STAGE_1[propellant_index]*1e6
    else:
        thrust_per_motor = THRUST_PER_MOTOR_STAGE_2[propellant_index]*1e6
        chamber_pressure = CHAMBER_PRESSURE_STAGE_2[propellant_index]*1e6
    n_motors = math.ceil(thrust_required/thrust_per_motor)
    return n_motors*237.8*(thrust_per_motor/chamber_pressure)**0.9375


def get_propellant_mass(stage, propellant_index, current_stage_mass, deltaV):
    # Leela
    # rocket equation, sea level isp for stage 1 and vacuum isp for stage 2
    if stage == 1:
        isp = ISP_SEA[propellant_index]
    else:
        isp = ISP_VAC[propellant_index]
    return current_stage_mass*(1 - math.exp(-deltaV/(isp*G0)))


# FUNCTION CATEGORY: checkers and post processing

def check_engines(n_engines, propellant_index, stage_index):
    # Colby
    # makes sure the engine nozzles physically fit under the stage. uses analyitical packing values for small counts, then switches to an area approximation
    exhaust_d = EXHAUST_DIAMETER[stage_index][propellant_index]
    stage_d = DIAMETER[stage_index]
    if n_engines in ENGINE_PACKING_FACTOR:
        if exhaust_d*ENGINE_PACKING_FACTOR[n_engines] > stage_d:
            raise RuntimeError("engine diameter exceeds limit")
    elif n_engines > 13:
        if n_engines*(exhaust_d/2)**2*math.pi > (stage_d/2)**2*math.pi*0.8:
            raise RuntimeError("engine diameter exceeds limit")


def print_vehicle(diameter1, diameter2, result):
    # Colby
    # prints the stage 1 / stage 2 mass breakdown side by side
    stage1_mass, stage2_mass, total_mass, cost, stage1_vars, stage2_vars = result
    print("diameter_1:", diameter1, "diameter_2:", diameter2)
    for name in stage1_vars:
        print(name, round(stage1_vars[name], 1), round(stage2_vars[name], 1))
    print("stage masses:", round(stage1_mass), round(stage2_mass))
    print("total mass:", round(total_mass), "cost:", round(cost))


if __name__ == "__main__":
    # 0 = LCH4/LOX, 1 = LH2/LOX, 2 = RP1/LOX, 3 = Solid, 4 = UDMH/N2O4
    stage_1_idx = 1
    stage_2_idx = 3
    deltaV1 = 7550
    results = []
    # sweep the diameter and keep every config that passes the constraints
    # for deltaV1 in np.linspace(1000, 10000, 100):
    for DIAMETER[0] in np.linspace(5.2, 40, 100):
        # for DIAMETER[1] in np.linspace(5.2, 20, 100):
        DIAMETER[1] = DIAMETER[0]
        try:
            results.append([DIAMETER[0], DIAMETER[1], size_vehicle(deltaV1, stage_1_idx, stage_2_idx), deltaV1])
        except RuntimeError:
            #can also catch errors to see what failed but don't really care now that debugging is done
            pass

    # pick out the lowest mass and lowest cost designs
    if not results:
        print("No valid configurations found.")
    else:
        min_mass_d1, min_mass_d2, min_mass_r, min_mass_dv1 = min(results, key=lambda r: r[2][2])
        print("\n\n\nMinimum mass results:")
        print("Delta V1:", min_mass_dv1)
        print_vehicle(min_mass_d1, min_mass_d2, min_mass_r)
        min_cost_d1, min_cost_d2, min_cost_r, min_cost_dv1 = min(results, key=lambda r: r[2][3])
        print("\n\n\nMinimum cost results:")
        print("Delta V1:", min_cost_dv1)
        print_vehicle(min_cost_d1, min_cost_d2, min_cost_r)
