
static FUEL_MIXTURE_MASS_RATIO: [f64; 5] = [3.6, 6.03, 2.72, 0.0, 2.67];
static THRUST_PER_MOTOR_STAGE_1: [f64; 5] = [2.26, 1.86, 1.92, 4.5, 1.75];
static EXHAUST_PER_MOTOR_STAGE_2: [f64; 5] = [0.745, 0.099, 0.061, 2.94, 0.067];
static EXHAUST_DIAMETER_STAGE_1: [f64; 5] = [2.4, 2.4, 3.7, 6.6, 1.5];
static EXHAUST_DIAMETER_STAGE_2: [f64; 5] = [1.5, 2.15, 0.92, 2.34, 1.13];
static CHAMBER_PRESSURE_STAGE_1: [f64; 5] = [35.16, 20.64, 25.8, 10.5, 15.7];
static CHAMBER_PRESSURE_STAGE_2: [f64; 5] = [10.1, 4.2, 6.77, 5.0, 14.7];
static NOZZLE_EXPANSION_STAGE_1: [f64; 5] = [34.34, 78.0, 37.0, 16.0, 26.2];
static NOZZLE_EXPANSION_STAGE_2: [f64; 5] = [45.0, 84.0, 14.5, 56.0, 81.3];
static RHO_PROPELLANT: [i64; 7] = [71, 1140, 820, 423, 1680, 1442, 791];
static SIGMA: [f64; 5] = [0.067, 0.075, 0.063, 0.087, 0.061];
static ISP_SEA: [f64; 5] = [327.0, 366.0, 311.0, 269.0, 285.0];
static ISP_VAC: [f64; 5] = [380.0, 452.0, 337.0, 279.0, 316.0];
static DELTAV_TOTAL: f64 = 12300.0;
static MPL: f64 = 26000.0;

fn main() {
    println!("Hello, world!");


}

fn 