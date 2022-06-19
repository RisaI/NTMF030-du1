#![feature(trait_alias)]

trait PotFunc = Fn(f64) -> f64;

fn ricatti_bessel(l: usize, x: f64) -> f64 {
    match l {
        0 => x.sin(),
        1 => x.sin() / x - x.cos(),
        2 => (3. / x.powi(2) - 1.) * x.sin() - 3. * x.cos() / x,
        3 => (x.powi(2) - 15.) * x.cos() / x.powi(2) + (15. - 6. * x.powi(2)) * x.sin() / x.powi(3),
        _ => unimplemented!(),
    }
}

fn ricatti_neumann(l: usize, x: f64) -> f64 {
    match l {
        0 => x.cos(),
        1 => x.cos() / x + x.sin(),
        2 => (3. / x.powi(2) - 1.) * x.cos() + 3. * x.sin() / x,
        3 => {
            -(x.powi(2) - 15.) * x.sin() / x.powi(2) + (15. - 6. * x.powi(2)) * x.cos() / x.powi(3)
        }
        _ => unimplemented!(),
    }
}

#[derive(Clone)]
struct SimParams {
    pub m: f64,
    pub l: usize,
    pub e: f64,

    pub steps: usize,
    pub r_max: f64,
}

impl SimParams {
    pub fn new(m: f64, l: usize, e: f64, steps: usize, r_max: f64) -> Self {
        Self {
            m,
            l,
            e,
            steps,
            r_max,
        }
    }

    #[allow(dead_code)]
    pub fn with_e(&self, e: f64) -> Self {
        let mut ret = self.clone();
        ret.e = e;
        ret
    }

    fn k_sqr(&self, r: f64, pot: &dyn PotFunc) -> f64 {
        if r == 0. {
            return 0.;
        }

        let l_sqr = (self.l * (1 + self.l)) as f64;

        2. * self.m * (self.e - pot(r) - l_sqr / (2. * self.m * r.powi(2)))
    }

    pub fn get_h(&self) -> f64 {
        self.r_max / (self.steps - 1) as f64
    }

    pub fn get_phase_shift(&self, result: &SimResult) -> f64 {
        let g = result.get_g();
        let k = (2. * self.m * self.e).sqrt();
        let l = self.l;

        let r_1 = self.r_max;
        let r_2 = self.r_max - self.get_h();

        ((ricatti_bessel(l, k * r_1) - g * ricatti_bessel(l, k * r_2))
            / (g * ricatti_neumann(l, k * r_2) - ricatti_neumann(l, k * r_1)))
        .atan()
    }

    pub fn get_partial_cross_section(&self, result: &SimResult) -> f64 {
        use std::f64::consts::PI;
        let k_sqr = 2. * self.m * self.e;

        4. * PI * (2. * self.l as f64 + 1.) * self.get_phase_shift(result).sin().powi(2) / k_sqr
    }
}

struct SimResult {
    pub wavefunc: Vec<f64>,
}

impl SimResult {
    pub fn get_g(&self) -> f64 {
        let w = &self.wavefunc;
        w[w.len() - 1] / w[w.len() - 2]
    }
}

fn calculate_wf<F: PotFunc>(params: &SimParams, pot: F) -> SimResult {
    let h = params.get_h();
    let mut pts = Vec::with_capacity(params.steps);

    pts.push(0.);
    pts.push(h); // A = h

    for i in 2..params.steps {
        let mut next = 2.
            * pts[i - 1]
            * (1. - (5. * h.powi(2) / 12.) * params.k_sqr((i - 1) as f64 * h, &pot));

        next -= pts[i - 2] * (1. + (h.powi(2) / 12.) * params.k_sqr((i - 2) as f64 * h, &pot));

        next /= 1. + (h.powi(2) / 12.) * params.k_sqr(i as f64 * h, &pot);

        // println!("{next}");
        pts.push(next);
    }

    SimResult { wavefunc: pts }
}

struct MatchinResult {
    pub left: Vec<f64>,
    pub right: Vec<f64>,
    pub h: f64,
}

impl MatchinResult {
    pub fn get_f(&self) -> f64 {
        (self.right[self.right.len() - 1] - self.left[self.left.len() - 2]) / self.h
    }

    pub fn into_wf(mut self) -> Vec<f64> {
        self.right.pop();
        self.right.pop();
        self.right.reverse();
        self.left.append(&mut self.right);

        self.left
    }
}

fn matching_func<F: PotFunc>(params: &SimParams, pot: F) -> MatchinResult {
    let h = params.get_h();
    let mut left = vec![0., h];

    for i in 2.. {
        let mut next = 2.
            * left[i - 1]
            * (1. - (5. * h.powi(2) / 12.) * params.k_sqr((i - 1) as f64 * h, &pot));

        next -= left[i - 2] * (1. + (h.powi(2) / 12.) * params.k_sqr((i - 2) as f64 * h, &pot));

        next /= 1. + (h.powi(2) / 12.) * params.k_sqr(i as f64 * h, &pot);

        if next < left[i - 1] || (i as f64 * h) > 0.95 {
            // Matching point
            break;
        } else {
            left.push(next);
        }
    }

    // let r_m = (left.len() - 1) as f64 * h; // Matching point r
    let right_count = (32. / h) as usize;
    let r_max = h * (left.len() + right_count - 2) as f64;
    let mut right = vec![0., h];

    for i in 2..=right_count {
        let mut next = 2.
            * right[i - 1]
            * (1. - (5. * h.powi(2) / 12.) * params.k_sqr(r_max - (i - 1) as f64 * h, &pot));

        next -= right[i - 2]
            * (1. + (h.powi(2) / 12.) * params.k_sqr(r_max - (i - 2) as f64 * h, &pot));

        next /= 1. + (h.powi(2) / 12.) * params.k_sqr(r_max - i as f64 * h, &pot);

        right.push(next);
    }

    let norm = left[left.len() - 1] / right[right.len() - 2];
    right.iter_mut().for_each(|v| {
        *v *= norm;
    });

    // if norm < 0. {
    // right.pop();
    // right.pop();
    // right.reverse();
    // left.append(&mut right);

    //     plot(
    //         "AAAA",
    //         [("asd".into(), (0..left.len()).map(|i| i as f64), left)],
    //     );

    //     panic!();
    // }

    MatchinResult { h, left, right }
}

fn task1_pot(a: f64, m: f64) -> impl PotFunc {
    move |r| {
        if r > a {
            0.
        } else {
            -1. * (4.8f64).powi(2) / (2. * m * a.powi(2))
        }
    }
}

fn task2_pot() -> impl PotFunc {
    use std::ops::Neg;

    move |r: f64| -3. * (r.powi(2) / 4.).neg().exp() + (r - 3.).powi(2).neg().exp()
}

fn plot(
    title: &str,
    traces: impl IntoIterator<
        Item = (
            String,
            impl IntoIterator<Item = f64>,
            impl IntoIterator<Item = f64>,
        ),
    >,
) {
    use plotly::{common::Title, Layout, Plot, Scatter};

    let mut plot = Plot::new();

    plot.set_layout(Layout::new().title(Title::new(title)));

    for (t_title, xs, ys) in traces.into_iter() {
        plot.add_trace(Scatter::new(xs, ys).name(t_title.as_str()));
    }

    plot.show();
}

fn find_root<F: PotFunc>(mut from: f64, mut to: f64, f: F) -> Option<f64> {
    while (to - from) > 1e-4 {
        let [a, b] = [from, to].map(&f);

        let midpoint = from - (to - from) * a / (b - a);

        if midpoint > to || midpoint < from {
            return None;
        }

        let c = f(midpoint);

        if (a * c).signum() < 0. {
            to = midpoint;
        } else {
            from = midpoint
        }
    }

    Some((from + to) / 2.)
}

fn equidist(from: f64, to: f64, elems: usize) -> impl Iterator<Item = (f64, f64)> + Clone {
    let d = (to - from) / elems as f64;
    (0..elems).map(move |i| (from + i as f64 * d, from + (i + 1) as f64 * d))
}

fn main() {
    let plots = true;

    // Task 1
    // Subtask 1 and 2 - above
    // Subtask 3
    if plots {
        let energies = (0..=100).map(|i| i as f64 * 5. * 1e-2);

        plot(
            "Zero potential",
            (0..=3).map(|l| {
                (
                    format!("$\\delta_{l}(E)$"),
                    energies.clone(),
                    energies.clone().map(move |e| {
                        let params = SimParams::new(1., l, e, 1_000, 10.);
                        let result = calculate_wf(&params, |_| 0.);

                        // println!("{}", result.get_g());

                        params.get_phase_shift(&result)
                    }),
                )
            }),
        );
    }

    // Subtask 4
    if plots {
        let energies = (0..=1000).map(|i| i as f64 * 50. * 1e-3);
        let results = (0..=3)
            .map(|l| {
                energies
                    .clone()
                    .map(|e| {
                        let params = SimParams::new(1., l, e, 10_000, 2.);
                        let result = calculate_wf(&params, &task1_pot(1., 1.));

                        // println!("{}", result.get_g());

                        (
                            params.get_phase_shift(&result),
                            params.get_partial_cross_section(&result),
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        plot(
            "Task 1 - phase shifts",
            results.iter().enumerate().map(|(l, r)| {
                (
                    format!("$\\delta_{l}(E)$"),
                    energies.clone(),
                    r.iter().map(|(a, _)| *a),
                )
            }),
        );

        plot(
            "Task 1 - cross sections",
            results.iter().enumerate().map(|(l, r)| {
                (
                    format!("$\\sigma_{l}(E)$"),
                    energies.clone(),
                    r.iter().map(|(_, b)| *b),
                )
            }),
        );
    }

    // Subtask 8 and 9
    if plots {
        plot(
            "Task 1 ",
            (0..=3).map(|l| {
                let params = SimParams {
                    l,
                    e: 0.,
                    m: 1.,
                    r_max: 1.,
                    steps: 1_000,
                };

                let min_energy =
                    task1_pot(1., 1.)(0.) + (l * (l + 1)) as f64 / (2. * params.m) + 0.1;
                println!("E_{{min}} = {min_energy} for l = {l}");
                let energies = (0..=1000).map(move |i| min_energy - min_energy * i as f64 * 1e-3);

                let mf = {
                    let pot = task1_pot(1., 1.);
                    move |e| matching_func(&params.with_e(e), &pot).get_f()
                };
                equidist(min_energy, 0., 4)
                    .filter_map(|(from, to)| find_root(from, to, &mf))
                    .enumerate()
                    .for_each(|(i, root)| {
                        println!("root = {root}");
                        let params = SimParams::new(1., l, root, 2_000, 16.);
                        let d = params.get_h();

                        plot(
                            &format!("Task 1, l = {l}"),
                            [(
                                format!("Bound state {}", i + 1),
                                (0..=params.steps).map(|i| i as f64 * d),
                                matching_func(&params, task1_pot(1., 1.)).into_wf(),
                            )],
                        );
                    });

                (format!("f_{l}(E)"), energies.clone(), energies.map(mf))
            }),
        );
    }

    // Task 2
    if plots {
        let energies = (0..=1000).map(|i| i as f64 * 1e-3);
        let results = (0..=1)
            .map(|l| {
                energies
                    .clone()
                    .map(|e| {
                        let params = SimParams::new(1., l, e, 10_000, 5.5);
                        let result = calculate_wf(&params, &task2_pot());

                        // println!("{}", result.get_g());

                        (
                            params.get_phase_shift(&result),
                            params.get_partial_cross_section(&result),
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        plot(
            "Task 2 - phase shifts",
            results.iter().enumerate().map(|(l, r)| {
                (
                    format!("$\\delta_{l}(E)$"),
                    energies.clone(),
                    r.iter().map(|(a, _)| *a),
                )
            }),
        );

        plot(
            "Task 2 - cross sections",
            results.iter().enumerate().map(|(l, r)| {
                (
                    format!("$\\sigma_{l}(E)$"),
                    energies.clone(),
                    r.iter().map(|(_, b)| *b),
                )
            }),
        );

        println!("Task 2 roots:");
        for l in 0..=1 {
            println!("l = {l}");
            let mf = {
                let pot = task2_pot();
                move |e| matching_func(&SimParams::new(1., l, e, 10_000, 5.5), &pot).get_f()
            };
            equidist(-3., 0., 4)
                .filter_map(|(from, to)| find_root(from, to, &mf))
                .enumerate()
                .for_each(|(_, root)| {
                    println!("root = {root}");
                });
        }
    }
}
