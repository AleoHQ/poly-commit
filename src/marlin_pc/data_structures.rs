use crate::{PCCommitment, PCCommitterKey, PCRandomness, PCVerifierKey, Vec};
use core::ops::{Add, AddAssign};
use rand_core::RngCore;
use snarkos_models::curves::PairingEngine;
use snarkos_utilities::bytes::{FromBytes, ToBytes};

use crate::kzg10;
/// `UniversalParams` are the universal parameters for the KZG10 scheme.
pub type UniversalParams<E> = kzg10::UniversalParams<E>;

/// `CommitterKey` is used to commit to and create evaluation proofs for a given
/// polynomial.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct CommitterKey<E: PairingEngine> {
    /// The key used to commit to polynomials.
    pub powers: Vec<E::G1Affine>,

    /// The key used to commit to shifted polynomials.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub shifted_powers: Option<Vec<E::G1Affine>>,

    /// The key used to commit to hiding polynomials.
    pub powers_of_gamma_g: Vec<E::G1Affine>,

    /// The degree bounds that are supported by `self`.
    /// In ascending order from smallest to largest.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub enforced_degree_bounds: Option<Vec<usize>>,
    /// The maximum degree supported by the `UniversalParams` `self` was derived
    /// from.
    pub max_degree: usize,
}

impl<E: PairingEngine> CommitterKey<E> {
    /// Obtain powers for the underlying KZG10 construction
    pub fn powers<'a>(&'a self) -> kzg10::Powers<'a, E> {
        kzg10::Powers {
            powers_of_g: self.powers.as_slice().into(),
            powers_of_gamma_g: self.powers_of_gamma_g.as_slice().into(),
        }
    }

    /// Obtain powers for committing to shifted polynomials.
    pub fn shifted_powers<'a>(
        &'a self,
        degree_bound: impl Into<Option<usize>>,
    ) -> Option<kzg10::Powers<'a, E>> {
        self.shifted_powers.as_ref().map(|shifted_powers| {
            let powers_range = if let Some(degree_bound) = degree_bound.into() {
                assert!(self
                    .enforced_degree_bounds
                    .as_ref()
                    .unwrap()
                    .contains(&degree_bound));
                let max_bound = self
                    .enforced_degree_bounds
                    .as_ref()
                    .unwrap()
                    .last()
                    .unwrap();
                (max_bound - degree_bound)..
            } else {
                0..
            };
            let ck = kzg10::Powers {
                powers_of_g: (&shifted_powers[powers_range]).into(),
                powers_of_gamma_g: self.powers_of_gamma_g.as_slice().into(),
            };
            ck
        })
    }
}

impl<E: PairingEngine> PCCommitterKey for CommitterKey<E> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.powers.len()
    }
}

impl<E: PairingEngine> FromBytes for CommitterKey<E> {
    fn read<R: snarkos_utilities::io::Read>(mut reader: R) -> snarkos_utilities::io::Result<Self> {
        let powers = Vec::<E::G1Affine>::read(&mut reader)?;
        let shifted_powers = Option::<Vec<E::G1Affine>>::read(&mut reader)?;
        let powers_of_gamma_g = Vec::<E::G1Affine>::read(&mut reader)?;
        let enforced_degree_bounds = match Option::<Vec<u64>>::read(&mut reader)? {
            Some(v) => Some(v.iter().map(|x| *x as usize).collect()),
            _ => None,
        };
        let max_degree = u64::read(&mut reader)? as usize;
        Ok(Self {
            powers,
            shifted_powers,
            powers_of_gamma_g,
            enforced_degree_bounds,
            max_degree,
        })
    }
}

impl<E: PairingEngine> ToBytes for CommitterKey<E> {
    fn write<W: snarkos_utilities::io::Write>(
        &self,
        mut writer: W,
    ) -> snarkos_utilities::io::Result<()> {
        self.powers.write(&mut writer)?;
        self.shifted_powers.write(&mut writer)?;
        self.powers_of_gamma_g.write(&mut writer)?;
        let transformed_enforced_degree_bounds = match &self.enforced_degree_bounds {
            Some(v) => Some(v.iter().map(|x| *x as u64).collect::<Vec<_>>()),
            _ => None,
        };
        transformed_enforced_degree_bounds.write(&mut writer)?;
        (self.max_degree as u64).write(&mut writer)
    }
}

/// `VerifierKey` is used to check evaluation proofs for a given commitment.
#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct VerifierKey<E: PairingEngine> {
    /// The verification key for the underlying KZG10 scheme.
    pub vk: kzg10::VerifierKey<E>,
    /// Information required to enforce degree bounds. Each pair
    /// is of the form `(degree_bound, shifting_advice)`.
    /// The vector is sorted in ascending order of `degree_bound`.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub degree_bounds_and_shift_powers: Option<Vec<(usize, E::G1Affine)>>,
    /// The maximum degree supported by the `UniversalParams` `self` was derived
    /// from.
    pub max_degree: usize,
    /// The maximum degree supported by the trimmed parameters that `self` is
    /// a part of.
    pub supported_degree: usize,
}

impl<E: PairingEngine> FromBytes for VerifierKey<E> {
    fn read<R: snarkos_utilities::io::Read>(mut reader: R) -> snarkos_utilities::io::Result<Self> {
        let vk = kzg10::VerifierKey::<E>::read(&mut reader)?;
        let is_exist = bool::read(&mut reader)?;
        let degree_bounds_and_shift_powers = if is_exist {
            let len = u64::read(&mut reader)?;
            let v = (0..len)
                .map(|_| {
                    let d = u64::read(&mut reader)? as usize;
                    let power = E::G1Affine::read(&mut reader)?;
                    Ok((d, power))
                })
                .collect::<snarkos_utilities::io::Result<Vec<_>>>()?;

            Some(v)
        } else {
            None
        };
        let max_degree = u64::read(&mut reader)? as usize;
        let supported_degree = u64::read(&mut reader)? as usize;

        Ok(Self {
            vk,
            degree_bounds_and_shift_powers,
            max_degree,
            supported_degree,
        })
    }
}

impl<E: PairingEngine> ToBytes for VerifierKey<E> {
    fn write<W: snarkos_utilities::io::Write>(
        &self,
        mut writer: W,
    ) -> snarkos_utilities::io::Result<()> {
        self.vk.write(&mut writer)?;
        match &self.degree_bounds_and_shift_powers {
            Some(v) => {
                true.write(&mut writer)?;
                (v.len() as u64).write(&mut writer)?;
                for (d, power) in v {
                    (*d as u64).write(&mut writer)?;
                    power.write(&mut writer)?;
                }
            }
            _ => {
                false.write(&mut writer)?;
            }
        }
        (self.max_degree as u64).write(&mut writer)?;
        (self.supported_degree as u64).write(&mut writer)
    }
}

impl<E: PairingEngine> VerifierKey<E> {
    /// Find the appropriate shift for the degree bound.
    pub fn get_shift_power(&self, bound: usize) -> Option<E::G1Affine> {
        self.degree_bounds_and_shift_powers.as_ref().and_then(|v| {
            v.binary_search_by(|(d, _)| d.cmp(&bound))
                .ok()
                .map(|i| v[i].1)
        })
    }
}

impl<E: PairingEngine> PCVerifierKey for VerifierKey<E> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.supported_degree
    }
}

/// Commitment to a polynomial that optionally enforces a degree bound.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Copy(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Commitment<E: PairingEngine> {
    pub(crate) comm: kzg10::Commitment<E>,
    pub(crate) shifted_comm: Option<kzg10::Commitment<E>>,
}

impl<E: PairingEngine> ToBytes for Commitment<E> {
    #[inline]
    fn write<W: snarkos_utilities::io::Write>(
        &self,
        mut writer: W,
    ) -> snarkos_utilities::io::Result<()> {
        self.comm.write(&mut writer)?;
        self.shifted_comm.write(&mut writer)
    }
}

impl<E: PairingEngine> FromBytes for Commitment<E> {
    #[inline]
    fn read<R: snarkos_utilities::io::Read>(mut reader: R) -> snarkos_utilities::io::Result<Self> {
        let comm = kzg10::Commitment::read(&mut reader)?;
        let shifted_comm = Option::<kzg10::Commitment<E>>::read(&mut reader)?;
        Ok(Self { comm, shifted_comm })
    }
}

impl<E: PairingEngine> PCCommitment for Commitment<E> {
    #[inline]
    fn empty() -> Self {
        Self {
            comm: kzg10::Commitment::empty(),
            shifted_comm: Some(kzg10::Commitment::empty()),
        }
    }

    fn has_degree_bound(&self) -> bool {
        self.shifted_comm.is_some()
    }

    fn size_in_bytes(&self) -> usize {
        self.comm.size_in_bytes() + self.shifted_comm.as_ref().map_or(0, |c| c.size_in_bytes())
    }
}

/// `Randomness` hides the polynomial inside a commitment. It is output by `KZG10::commit`.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Randomness<E: PairingEngine> {
    pub(crate) rand: kzg10::Randomness<E>,
    pub(crate) shifted_rand: Option<kzg10::Randomness<E>>,
}

impl<'a, E: PairingEngine> Add<&'a Self> for Randomness<E> {
    type Output = Self;

    fn add(mut self, other: &'a Self) -> Self {
        self += other;
        self
    }
}

impl<'a, E: PairingEngine> AddAssign<&'a Self> for Randomness<E> {
    #[inline]
    fn add_assign(&mut self, other: &'a Self) {
        self.rand += &other.rand;
        if let Some(r1) = &mut self.shifted_rand {
            *r1 += other
                .shifted_rand
                .as_ref()
                .unwrap_or(&kzg10::Randomness::empty());
        } else {
            self.shifted_rand = other.shifted_rand.as_ref().map(|r| r.clone());
        }
    }
}

impl<'a, E: PairingEngine> Add<(E::Fr, &'a Randomness<E>)> for Randomness<E> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: (E::Fr, &'a Randomness<E>)) -> Self {
        self += other;
        self
    }
}

impl<'a, E: PairingEngine> AddAssign<(E::Fr, &'a Randomness<E>)> for Randomness<E> {
    #[inline]
    fn add_assign(&mut self, (f, other): (E::Fr, &'a Randomness<E>)) {
        self.rand += (f, &other.rand);
        let empty = kzg10::Randomness::empty();
        if let Some(r1) = &mut self.shifted_rand {
            *r1 += (f, other.shifted_rand.as_ref().unwrap_or(&empty));
        } else {
            self.shifted_rand = other.shifted_rand.as_ref().map(|r| empty + (f, r));
        }
    }
}

impl<E: PairingEngine> PCRandomness for Randomness<E> {
    fn empty() -> Self {
        Self {
            rand: kzg10::Randomness::empty(),
            shifted_rand: None,
        }
    }

    fn rand<R: RngCore>(hiding_bound: usize, has_degree_bound: bool, rng: &mut R) -> Self {
        let shifted_rand = if has_degree_bound {
            Some(kzg10::Randomness::rand(hiding_bound, false, rng))
        } else {
            None
        };
        Self {
            rand: kzg10::Randomness::rand(hiding_bound, false, rng),
            shifted_rand,
        }
    }
}

impl<E: PairingEngine> FromBytes for Randomness<E> {
    fn read<R: snarkos_utilities::io::Read>(mut reader: R) -> snarkos_utilities::io::Result<Self> {
        let rand = kzg10::Randomness::<E>::read(&mut reader)?;
        let shifted_rand = Option::<kzg10::Randomness<E>>::read(&mut reader)?;
        Ok(Self { rand, shifted_rand })
    }
}

impl<E: PairingEngine> ToBytes for Randomness<E> {
    fn write<W: snarkos_utilities::io::Write>(
        &self,
        mut writer: W,
    ) -> snarkos_utilities::io::Result<()> {
        self.rand.write(&mut writer)?;
        self.shifted_rand.write(&mut writer)
    }
}
