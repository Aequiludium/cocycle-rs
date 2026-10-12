//! Read owned F3 persistence data and representatives after releasing the source.
use cocycle::algebra::PrimeField;
use cocycle::descriptors::{betti_curve, finite_lifetime_summary};
use cocycle::diagram::{PersistenceData, Representative, RepresentativeKind};
use cocycle::diagram_distances::bottleneck_distance_results;
use cocycle::filtration::RipsBuilder;
use cocycle::geometry::{DissimilarityMatrixView, MatrixLayout};
use cocycle::persistence::{PersistenceExt, RepresentativeRequest, RepresentativeSelection};

fn main() -> cocycle::Result<()> {
    let result = {
        // Four-cycle edges enter at 1; diagonals enter at 2 and fill the cycle.
        let values = [1., 2., 1., 1., 2., 1.];
        let matrix = DissimilarityMatrixView::new(&values, 4, MatrixLayout::LowerTriangle)?;
        let requests = [RepresentativeRequest::new(
            1,
            1.,
            RepresentativeSelection::Both,
        )?];
        RipsBuilder::from_distance_matrix(matrix)
            .persistence()
            .field(PrimeField::new(3)?)
            .representatives(&requests)
            .compute()?
    }; // Input buffer, builder and requests are no longer available.
    let borrowed: &PersistenceData = result.as_ref();
    assert!(std::ptr::eq(borrowed.diagram(), result.diagram()));
    let (data, representatives) = result.into_parts(); // Move; keep witness ordinals.
    assert_eq!(betti_curve(data.diagram(), 1, &[1., 2.])?, [1, 0]);
    assert_eq!(
        finite_lifetime_summary(data.diagram(), 1)?.total_persistence(),
        1.
    );
    assert_eq!(bottleneck_distance_results(&data, &data, 1)?, 0.);
    read_representatives(&data, representatives.as_deref().unwrap());
    Ok(())
}

// This consumer needs only interpreted results, with no builder or reducer access.
fn read_representatives(data: &PersistenceData, representatives: &[Representative]) {
    assert_eq!(representatives.len(), 2);
    assert_eq!(data.context().characteristic(), 3);
    for representative in representatives {
        let interval = data
            .diagram()
            .interval(representative.interval_index())
            .unwrap();
        let kind = match representative.kind() {
            RepresentativeKind::Cycle => "cycle",
            RepresentativeKind::Cocycle => "cocycle",
        };
        println!(
            "{kind} at {} for interval {:?} over F{}:",
            representative.scale(),
            interval,
            representative.characteristic()
        );
        for term in representative.terms() {
            println!("  {:?}: {}", term.vertices(), term.coefficient());
        }
    }
}
