use peroxide::fuga::*;
use std::collections::HashMap;

fn main() {
    let m = 1f64;
    let k = 200f64;
    let c_s = 2f64 * (k * m).sqrt();
    let zeta_vec = vec![0f64, 0.01f64, 0.02f64];
    let c_vec = zeta_vec.fmap(|zeta| c_s * zeta);

    let x_init = 0.1f64;
    let v_init = 0f64;
    let a_init = -20f64;

    let dt = 1e-3;

    let gamma = 0.5;
    let beta = 0.25; // for average constant acceleration
    
    let mut damped_sho_vec = c_vec.into_iter()
        .map(|c| NewmarkSHO::new(x_init, v_init, a_init, m, c, k, gamma, beta))
        .collect::<Vec<NewmarkSHO>>();

    let t_vec = seq(0f64, 10f64, dt);
    let mut x_hash = HashMap::new();
    let mut v_hash = HashMap::new();
    let mut a_hash = HashMap::new();

    for i in 0 .. zeta_vec.len() {
        x_hash.insert(i, vec![0f64; t_vec.len()]);
        v_hash.insert(i, vec![0f64; t_vec.len()]);
        a_hash.insert(i, vec![0f64; t_vec.len()]);
    }

    for i in 0 .. t_vec.len() {
        for j in 0 .. zeta_vec.len() {
            let (x, v, a) = damped_sho_vec[j].get_state();
            let x_vec = x_hash.get_mut(&j).unwrap();
            let v_vec = v_hash.get_mut(&j).unwrap();
            let a_vec = a_hash.get_mut(&j).unwrap();

            x_vec[i] = x;
            v_vec[i] = v;
            a_vec[i] = a;

            damped_sho_vec[j].step(dt);
        }
    }

    let mut df = DataFrame::new(vec![]);
    df.push("t", Series::new(t_vec));
    for i in 0 .. zeta_vec.len() {
        df.push(&format!("x_{}", i), Series::new(x_hash[&i].clone()));
        df.push(&format!("v_{}", i), Series::new(v_hash[&i].clone()));
        df.push(&format!("a_{}", i), Series::new(a_hash[&i].clone()));
        df.push(&format!("zeta_{}", i), Series::new(vec![zeta_vec[i]]));
    }

    df.print();

    df.write_nc("data/sho_newmark_beta.nc").expect("Can't write newmark_beta");
}

pub struct NewmarkSHO {
    x: f64,     // position
    v: f64,     // velocity
    a: f64,     // acceleration
    m: f64,     // mass
    c: f64,     // damping coefficient
    k: f64,     // elastic constant
    gamma: f64, // newmark gamma
    beta: f64,  // newmark beta
}

impl NewmarkSHO {
    pub fn new(x: f64, v: f64, a: f64, m: f64, c: f64, k: f64, gamma: f64, beta: f64) -> NewmarkSHO {
        NewmarkSHO {
            x,
            v,
            a,
            m,
            c,
            k,
            gamma,
            beta,
        }
    }

    pub fn get_state(&self) -> (f64, f64, f64) {
        (self.x, self.v, self.a)
    }

    pub fn step(&mut self, dt: f64) {
        let x = self.x;
        let v = self.v;
        let a = self.a;
        let m = self.m;
        let c = self.c;
        let k = self.k;
        let gamma = self.gamma;
        let beta = self.beta;

        let a_next = - (k*x + (c + k*dt)*v + ((1f64-gamma)*c*dt + 0.5 * (1f64-2f64*beta)*k*dt.powi(2))*a) / (m + c*gamma*dt + k*beta*dt.powi(2));
        let v_next = v + (1f64-gamma)*dt*a + gamma*dt*a_next;
        let x_next = x + dt*v + 0.5*dt.powi(2)*((1f64-2f64*beta)*a + 2f64*beta*a_next);

        self.x = x_next;
        self.v = v_next;
        self.a = a_next;
    }
}
