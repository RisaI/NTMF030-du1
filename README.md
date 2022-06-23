# Homework 1

## Task 1
### Subtask 1, 2 and 3
The implementation can be seen in `src/main.rs`.

### Subtask 4 and 5
The following phase shifts and partial cross-sections were calculated:
![](img/t1_phase-shift.png)
![](img/t1_cross-sections.png)

The Ramsauer-Townsend effect can be observed around $E \approx 3.6$ for $l = 0$, where the phase shift crosses zero and the potential becomes opaque for the scattering particle, not affecting it in a measurable way. The cross-section approaches zero in this region.

A sharp resonance can be seen around $E \approx 3.3$ for $l = 3$. The centrifugal term enables a state with positive energy which would be bound in classical mechanics, but since QM permits quantum tunneling, the state is decaying. This can be detected both as a delay between projectile release and detection and a sharp peak in the cross-section diagram.


### Subtask 6 and 7
Interesting case was the case of $ l = 2,3 $. For $l = 2$, a resonance could be achieved by making the potential slightly weaker (by a factor around 0.66). For $l = 3$, the already existing resonance can be destroyed by tuning the potential.

![](img/t1_l2_ps.png)
![](img/t1_l3_ps.png)

Since the centrifugal term is zero, no resonance could be achieved by tuning the potential strength for $l = 0$:
![](img/t1_l0_ps.png)
![](img/t1_l0_cs.png)

### Subtask 8 and 9
The matching functions are the following:
![](img/t1_f0.png)
![](img/t1_f1.png)
![](img/t1_f2.png)
![](img/t1_f3.png)

Roots were found using a bisection algorithm. For $l = 0$, there were two roots corresponding with bound states (-8.2066, -0.0666) and one false positive (-4.3117). For $l = 1$, there was one bound state (-4.88168) and one false positive (-0.31279), and for $l = 2$ there was a single root (-0.96017). $l = 3$ had no roots.

### Subtask 10 and 11
![](img/t1_bound-states.png)

The peaks of the bound states are moving into higher radial distance $r$ with increasing $l$. Their energy is also steadily increasing, and it can be inferred that the bound state for $l = 3$ would have a positive energy. This is precisely the resonance we have described earlier.

This directly ties to the _Levinson's theorem_ which states that $\delta_l(0) - \delta_l(\infty) = n_l \pi$, where $n_l$ is the number of bound states for $l$. It can be seen in the very first graph, that the phase shift for $l=3$ loops upward and then back again, then converging to zero in the same modulo segment as the initial point. $l=2$ loops downward twice, as it has two bound states and $l=1,2$ only once.

---


## Task 2
![](img/t2_phase-shift.png)

In the phase shift diagram, a resonance is visible for $l = 1$ and $E \approx 0.512$. This also emerges as a noticeable bump in the cross-section diagram:

![](img/t2_cross-sections.png)

Bound states were determined the same way as in task 1. There was one false-positive and two physical bound states.

![](img/t2_f.png)
![](img/t2_bound-states.png)
