use jubjub::*;
use pairing::*;

#[derive(Clone, Copy)]
pub enum Personalization {
    NoteCommitment,
    MerkleTree(usize)
}

impl Personalization {
    pub fn get_bits(&self) -> Vec<bool> {
        match *self {
            Personalization::NoteCommitment =>
                vec![true, true, true, true, true, true],
            Personalization::MerkleTree(num) => {
                assert!(num < 63);

                (0..6).map(|i| (num >> i) & 1 == 1).collect()
            }
        }
    }
}

pub fn pedersen_hash<E, I>(
    personalization: Personalization,
    bits: I,
    params: &E::Params
) -> edwards::Point<E, PrimeOrder>
    where I: IntoIterator<Item=bool>,
          E: JubjubEngine
{
    let mut bits = personalization.get_bits().into_iter().chain(bits.into_iter());

    let mut result = edwards::Point::zero();
    let mut generators = params.pedersen_hash_generators().iter();

    loop {
        let mut acc = E::Fs::zero();
        let mut cur = E::Fs::one();
        let mut chunks_remaining = params.pedersen_hash_chunks_per_generator();
        let mut encountered_bits = false;

        // Grab three bits from the input
        // Do this 63 times
        while let Some(a) = bits.next() {
            encountered_bits = true;

            // Pad length to a multiple of 3
            let b = bits.next().unwrap_or(false);
            let c = bits.next().unwrap_or(false);

            // Start computing this portion of the scalar
            // Comments show the extreme values during the last iteration.
            let mut tmp = cur;          // tmp = 2^248  (62 * 4 doublings)
            if a {
                tmp.add_assign(&cur);   // 2^248 <= tmp <= 2^249
            }
            cur.double();               // cur = 2^249
            if b {
                tmp.add_assign(&cur);   // tmp <= 2^250
            }

            // conditionally negate
            if c {
                tmp.negate();           // -2^250 <= tmp
            }

            // So far, -2^247 < acc < 2^247 because the largest tmp has been 2^246.

            // Accumulate this chunk.
            acc.add_assign(&tmp);

            // After the last iteration, and wrapping negative values modulo s:
            // 0  <  largest positive acc  <  2^250 + 2^247  <  s - 2^250 - 2^247)  <  lowest negative acc mod s  <  s

            chunks_remaining -= 1;

            if chunks_remaining == 0 {
                break;
            } else {
                cur.double(); // 2^2 * cur
                cur.double(); // 2^3 * cur
                cur.double(); // 2^4 * cur
            }
        }

        if !encountered_bits {
            break;
        }

        let mut tmp = generators.next().expect("we don't have enough generators").clone();
        tmp = tmp.mul(acc, params);
        result = result.add(&tmp, params);
    }

    result
}

#[cfg(test)]
pub mod test {

    use pairing::bls12_381::{Bls12};
    use super::*;

    pub struct TestVector<'a> {
        pub personalization: Personalization,
        pub input_bits: Vec<u8>,
        pub hash_x: &'a str,
        pub hash_y: &'a str,
    }

    use tests::pedersen_hash_vectors;

    #[test]
    fn test_pedersen_hash_points() {

        let test_vectors = pedersen_hash_vectors::get_vectors();

        let params = &JubjubBls12::new();

        let v = &test_vectors[0];
        let input_bools: Vec<bool> = v.input_bits.iter().map(|&i| i == 1).collect();

        // The 6 bits prefix is handled separately
        assert_eq!(v.personalization.get_bits(), &input_bools[..6]);

        let (x, y) = pedersen_hash::<Bls12, _>(
            v.personalization,
            input_bools.into_iter().skip(6),
            params,
        ).into_xy();

        assert_eq!(x.to_string(), v.hash_x);
        assert_eq!(y.to_string(), v.hash_y);
    }
}
