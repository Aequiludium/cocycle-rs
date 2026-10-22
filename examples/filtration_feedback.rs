//! Application-owned recomputation on one analytically differentiable region.
//! This is a finite experiment, not an autodiff or training API.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::descriptors::finite_lifetime_summary;
use cocycle::diagram::PersistenceDiagram;
use cocycle::persistence::PersistenceExt;
use cocycle::{Error, Result};

fn compute(theta: f64) -> Result<PersistenceDiagram> {
    if !theta.is_finite() {
        return Err(Error::NonFiniteValue {
            field: "feedback theta",
            index: None,
        });
    }
    let source = SimplicialComplex::new(vec![
        Simplex::new(vec![0], 0.)?,
        Simplex::new(vec![1], 0.)?,
        Simplex::new(vec![0, 1], theta.max(0.))?,
    ])?;
    Ok(source
        .persistence()
        .max_homology_dimension(0)
        .compute()?
        .diagram()
        .clone())
}

fn loss(theta: f64) -> Result<f64> {
    let mass = finite_lifetime_summary(&compute(theta)?, 0)?.total_persistence();
    Ok(0.5 * (mass - 1.).powi(2))
}

fn gradient(theta: f64) -> Result<f64> {
    if !theta.is_finite() || theta <= 0. {
        return Err(Error::InvalidParameter {
            parameter: "feedback gradient",
            reason: "this experiment differentiates only the region theta > 0",
        });
    }
    // On this region there is exactly one finite H0 bar [0,theta).
    Ok(theta - 1.)
}

fn main() -> Result<()> {
    let original = compute(2.)?;
    let mut theta = 2.;
    let h = 1. / 1024.;
    for iteration in 0..=8 {
        let diagram = compute(theta)?; // New source/result for this sample key.
        let summary = finite_lifetime_summary(&diagram, 0)?;
        let value = loss(theta)?;
        let derivative = gradient(theta)?;
        let difference = (loss(theta + h)? - loss(theta - h)?) / (2. * h);
        assert_eq!(derivative, difference); // All fixture values are dyadic.
        println!(
            "{{\"kind\":\"feedback\",\"snapshot\":\"edge-{iteration}\",\"iteration\":{iteration},\"theta\":{theta},\"mass\":{},\"finite_count\":{},\"loss\":{value},\"gradient\":{derivative},\"difference\":{difference},\"h\":{h}}}",
            summary.total_persistence(),
            summary.finite_count(),
        );
        theta -= 0.5 * derivative;
    }
    // Retained results keep the meaning of their old source after updates.
    assert_eq!(
        finite_lifetime_summary(&original, 0)?.total_persistence(),
        2.
    );
    for theta in [-h, 0., h] {
        let diagram = compute(theta)?;
        let summary = finite_lifetime_summary(&diagram, 0)?;
        println!(
            "{{\"kind\":\"tie\",\"theta\":{theta},\"mass\":{},\"finite_count\":{},\"loss\":{}}}",
            summary.total_persistence(),
            summary.finite_count(),
            loss(theta)?,
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recomputed_ph_loss_follows_the_analytic_contraction() -> Result<()> {
        let original = compute(2.)?;
        let mut theta = 2.;
        for step in 0..=8 {
            assert_eq!(theta, 1. + 0.5_f64.powi(step));
            assert_eq!(loss(theta)?, 0.5 * 0.25_f64.powi(step));
            for h in [1. / 16., 1. / 64., 1. / 1024.] {
                assert_eq!(
                    (loss(theta + h)? - loss(theta - h)?) / (2. * h),
                    gradient(theta)?
                );
            }
            theta -= 0.5 * gradient(theta)?;
        }
        assert_eq!(
            finite_lifetime_summary(&original, 0)?.total_persistence(),
            2.
        );
        Ok(())
    }

    #[test]
    fn zero_lifetime_tie_has_distinct_one_sided_slopes_and_no_smooth_gradient() -> Result<()> {
        for h in [1. / 16., 1. / 64., 1. / 1024.] {
            assert_eq!((loss(h)? - loss(0.)?) / h, -1. + h / 2.);
            assert_eq!((loss(0.)? - loss(-h)?) / h, 0.);
        }
        assert_eq!(finite_lifetime_summary(&compute(0.)?, 0)?.finite_count(), 0);
        assert_eq!(loss(-1.)?, 0.5); // Flat region is not the global minimum.
        for theta in [0., -1., f64::NAN, f64::INFINITY] {
            assert!(gradient(theta).is_err());
        }
        for theta in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(compute(theta).is_err());
        }
        Ok(())
    }
}
