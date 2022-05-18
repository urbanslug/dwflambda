/*!
Configs and related WFA types
 */

use super::utils;
use std::ops::Index;

// ---------
// Constants
// ---------

const NULL_OFFSET: i32 = -10;

// ----------------------
//         Config
// ----------------------
// TODO: use u8
pub struct Penalties {
    pub mismatch: i32,
    pub matches: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

pub struct Config {
    pub adapt: bool,
    // pub segment_length: u32, // segment size in bytes
    // pub step_size: u32,
    // pub thread_count: usize,
    pub verbosity: u8,
    pub penalties: Penalties,
}

// ----------------------
//     Types
// ----------------------
// WF Next
#[derive(Debug)]
pub struct AWFSet<'a> {
    // In
    pub in_m_sub: Option<&'a WaveFront>,
    pub in_m_gap: Option<&'a WaveFront>,
    pub in_i_ext: Option<&'a WaveFront>,
    pub in_d_ext: Option<&'a WaveFront>,

    // out
    pub out_m: Option<&'a WaveFront>,
    pub out_i: Option<&'a WaveFront>,
    pub out_d: Option<&'a WaveFront>,
}

// ----------------------
//     Core types
// ----------------------
#[derive(Debug, PartialEq, Eq)]
pub enum WfType {
    D,
    I,
    M,
}

#[derive(Debug, PartialEq, Eq)]
pub enum BacktraceOperation {
    MatchMismatch,
    Insertion,
    Deletion,
}

// TODO: should all be i32
// matrix offset, text offset & query offset
// pub type Offset = Vec<i32>;
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Offset {
    pub data: Vec<i32>,
    pub abdandoned: Option<Vec<u8>>,
}

impl Offset {
    // ---------
    // Construct
    // ---------
    pub fn blank() -> Self {
        Self {
            data: vec![0],
            abdandoned: None,
        }
    }

    // rename from_value
    pub fn from([e]: [i32; 1]) -> Self {
        Self {
            data: vec![e],
            abdandoned: None,
        }
    }

    pub fn from_vec(data: &Vec<i32>) -> Self {
        let z = data.len();
        Self {
            data: data.clone(),
            abdandoned: Some(vec![0; z]),
        }
    }

    pub fn null() -> Self {
        Self {
            data: vec![NULL_OFFSET],
            abdandoned: None,
        }
    }

    pub fn abandoned_all_null(&mut self) {
        let z = self.data.len();
        self.abdandoned = Some(vec![0; z]);
    }

    pub fn first(&self) -> i32 {
        self.data[0]
    }

    pub fn get_mut(&mut self, index: usize) -> &mut i32 {
        self.data.get_mut(index).expect("types::Offset::get_mut")
    }

    pub fn get_mut_maybe(&mut self, index: usize) -> Option<&mut i32> {
        self.data.get_mut(index)
    }

    // TODO: make this a result type
    pub fn set(&mut self, index: usize, value: i32) {
        self.data[index] = value
    }

    // TODO: make this a result type
    pub fn push(&mut self, value: i32) {
        self.data.push(value);
    }
    pub fn max(&self) -> i32 {
        *self.data.iter().max().unwrap()
    }

    pub fn min(&self) -> i32 {
        *self.data.iter().min().unwrap()
    }

    pub fn offset_count(&self) -> usize {
        self.data.len()
    }

    pub fn set_abandon(&mut self, index: usize) {
        self.abdandoned.as_mut().map(|v: &mut Vec<u8>| v[index] = 1);
    }

    pub fn is_abandoned(&self, index: usize) -> bool {
        match &self.abdandoned {
            None => false,
            Some(v) => {
                if v[index] == 1 {
                    true
                } else {
                    false
                }
            }
        }
    }
}

impl Index<usize> for Offset {
    type Output = i32;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.data[idx]
    }
}

/// The a single wavefront with a score
/// The furthest reaching point of a single wavefront
#[derive(Debug, Clone)]
pub struct WaveFront {
    /// the highest diagonal touched by the wavefront
    pub hi: i32,

    /// the lowest diagonal touched by the wavefront
    pub lo: i32,

    /// The offsets of each diagonal between hi and lo
    /// vals\[0\] is the score of the wavefront at diagonal hi and
    /// vals\[<last>\] is the score of the wavefront at diagonal lo
    /// length is (hi - lo) + 1
    pub offsets: Vec<Offset>,
}

impl WaveFront {
    pub fn new(hi: i32, lo: i32) -> Self {
        let len = utils::new_compute_wave_length(lo, hi);

        let starting_offset = Offset::from([0]);

        Self {
            hi,
            lo,
            offsets: vec![starting_offset; len],
        }
    }

    pub fn len(&self) -> usize {
        // TODO merge with utils

        let wave_len = utils::new_compute_wave_length(self.lo, self.hi);
        // TODO: remove
        assert_eq!(wave_len, self.offsets.len());

        wave_len
    }

    pub fn k_index(&self, k: i32) -> usize {
        utils::new_compute_k_index(k, self.lo, self.hi)
    }

    pub fn in_bounds(&self, k: i32) -> bool {
        utils::k_in_bounds(k, self.lo, self.hi)
    }

    // TODO: just return the i32
    /// return the offset at the k diagonal
    /// Computes the k-index internally
    /// takes the diagonal k (not k-index)
    pub fn get_offset(&self, k: i32) -> Option<&Offset> {
        if !utils::k_in_bounds(k, self.lo, self.hi) {
            // eprintln!("out of bounds");
            return None;
        }

        let k_index = utils::new_compute_k_index(k, self.lo, self.hi);

        self.offsets
            .get(k_index)
            .and_then(|offset: &Offset| Some(offset))
    }

    // TODO: rename to get_matrix offset
    pub fn get_offset_mut(&mut self, k: i32) -> Option<&mut Offset> {
        if !utils::k_in_bounds(k, self.lo, self.hi) {
            // eprintln!("out of bounds mut");
            return None;
        }

        let k_index = utils::new_compute_k_index(k, self.lo, self.hi);
        self.offsets
            .get_mut(k_index)
            .and_then(|offset: &mut Offset| Some(offset))
    }

    // Actual get offset
    pub fn bar(&self, k: i32) -> Option<&Offset> {
        if !utils::k_in_bounds(k, self.lo, self.hi) {
            // eprintln!("out of bounds mut");
            return None;
        }

        let k_index = utils::new_compute_k_index(k, self.lo, self.hi);
        self.offsets
            .get(k_index)
            .and_then(|offset: &Offset| Some(offset))
    }

    // Actual get offset mut
    pub fn foo(&mut self, k: i32) -> Option<&mut Offset> {
        if !utils::k_in_bounds(k, self.lo, self.hi) {
            // eprintln!("out of bounds mut");
            return None;
        }

        let k_index = utils::new_compute_k_index(k, self.lo, self.hi);
        self.offsets
            .get_mut(k_index)
            .and_then(|offset: &mut Offset| Some(offset))
    }
}

/// The set of wavefronts at a certain score
#[derive(Debug, Clone)]
pub struct WaveFrontSet {
    /// insertion wavefront
    pub i: Option<WaveFront>,

    /// deletion wavefront
    pub d: Option<WaveFront>,

    /// match wavefront
    /// $\tilde{M}_{s,k}$ is the value of the m wavefront at diagonal k
    pub m: Option<WaveFront>,
}

/// All the wavefronts
#[derive(Clone)]
pub struct WaveFronts {
    /// The set of wavefronts with each score, the index represents the score
    /// and, each element is a wavefront.
    /// WF_s is wavefront_set\[s\]
    pub wavefront_set: Vec<Option<WaveFrontSet>>,

    pub min_k: isize, // -qlen
    pub max_k: isize, // tlen
    pub a_k: i32,
}

impl WaveFronts {
    /// The scores should always be positive numbers
    pub fn get(&self, score: usize) -> &Option<WaveFrontSet> {
        &self.wavefront_set[score]
    }

    pub fn option_get(&self, score: usize) -> Option<&WaveFrontSet> {
        self.wavefront_set
            .get(score)
            .and_then(|maybe_wf_set| maybe_wf_set.as_ref())
    }

    pub fn get_wavefronts(&self, score: usize) -> Option<&WaveFrontSet> {
        self.wavefront_set
            .get(score)
            .and_then(|maybe_wf_set| maybe_wf_set.as_ref())
    }

    pub fn get_m_wavefront(&self, score: i32) -> Option<&WaveFront> {
        if score < 0 {
            return None;
        }

        let score = score as usize;
        let maybe_wf_set: Option<&WaveFrontSet> = self.option_get(score);
        match maybe_wf_set {
            Some(v) => v.m.as_ref(),
            _ => None,
        }
    }

    pub fn get_i_wavefront(&self, score: i32) -> Option<&WaveFront> {
        if score < 0 {
            return None;
        }

        let score = score as usize;
        let maybe_wf_set: Option<&WaveFrontSet> = self.option_get(score);
        match maybe_wf_set {
            Some(v) => v.i.as_ref(),
            _ => None,
        }
    }

    pub fn get_d_wavefront(&self, score: i32) -> Option<&WaveFront> {
        if score < 0 {
            return None;
        }

        let score = score as usize;
        let maybe_wf_set: Option<&WaveFrontSet> = self.option_get(score);
        match maybe_wf_set {
            Some(v) => v.d.as_ref(),
            _ => None,
        }
    }

    pub fn set_i_d_m(&mut self, score: usize) {
        let wf_set: &mut Option<WaveFrontSet> = self.wavefront_set.get_mut(score).unwrap();
        let wf_set: &mut Option<WaveFrontSet> = &mut self.wavefront_set[score];
    }

    pub fn len(&self) -> usize {
        self.wavefront_set.len()
    }

    pub fn max_score(&self) -> u32 {
        (self.len() - 1) as u32
    }

    // Allocate wavefronts (I, D, or M) for the given score
    pub fn allocate_wavefronts(
        &mut self,
        score: u32,
        lo: i32,
        hi: i32,
        wavefronts_to_allocate: &Vec<WfType>,
    ) -> Result<(), &str> {
        // should only add what is necessary
        let max_score = self.max_score();
        let len = num::abs_sub(hi, lo) as usize + 1;

        if max_score >= score {
            // we are trying to add a score that exists
            // eprintln!("previous score {} score {}", prev_score, score);
            return Err("[types::allocate_wavefronts] fuckery detected");
        }

        for index in max_score + 1..=score {
            if index == score {
                let wf_set = WaveFrontSet {
                    i: {
                        if wavefronts_to_allocate.contains(&WfType::I) {
                            Some(WaveFront::new(hi, lo))
                        } else {
                            None
                        }
                    },
                    d: {
                        if wavefronts_to_allocate.contains(&WfType::D) {
                            Some(WaveFront::new(hi, lo))
                        } else {
                            None
                        }
                    },
                    m: {
                        if wavefronts_to_allocate.contains(&WfType::M) {
                            Some(WaveFront::new(hi, lo))
                        } else {
                            None
                        }
                    },
                };

                self.wavefront_set.push(Some(wf_set));
            } else {
                self.wavefront_set.push(None);
            }
        }
        Ok(())
    }
}
