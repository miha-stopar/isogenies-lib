use thiserror::Error;

#[derive(Clone, Debug, Eq, Error, PartialEq)]
#[error("Not quadratic residue")]
pub struct NotQuadraticResidueError();

#[derive(Clone, Debug, Eq, Error, PartialEq)]
#[error("Not coprime")]
pub struct NotCoprime();

#[derive(Clone, Debug, Eq, Error, PartialEq)]
#[error("Equivalent prime ideal not found")]
pub struct EquivalentPrimeIdealNotFound();

#[derive(Clone, Debug, Eq, Error, PartialEq)]
#[error("Represent integer not found")]
pub struct RepresentIntegerNotFound();

#[derive(Clone, Debug, Eq, Error, PartialEq)]
#[error("Short dimension 2 vector not found")]
pub struct ShortDim2VecNotFound();

#[derive(Clone, Debug, Eq, Error, PartialEq)]
#[error("KLPT not successful")]
pub struct KLPTNotSuccessful();
