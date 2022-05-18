/*!
Abstract WFA using match and traceback functions
 */
use std::cmp;

use super::types;
use super::utils::{self};
use super::wf_extend::wf_extend;
use super::wf_next as core;
use super::wf_traceback::wf_traceback;
use fbox::macros::max;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

// TODO: return Result type
pub fn wf_align<F, G>(
    tlen: u32,
    qlen: u32,
    config: &types::Config,
    match_lambda: &mut F,
    traceback_lambda: &mut G,
) -> Result<(usize, String), String>
where
    F: FnMut(&mut i32, &mut i32, &mut types::Offset) -> bool,
    G: FnMut((i32, i32), (i32, i32)) -> bool,
{
    if config.verbosity > 1 {
        eprintln!("[wflambda::wf_align]");
    }

    // compute the central diagonal, a_k.
    let a_k: i32 = (tlen as isize - qlen as isize) as i32;

    // the furthest offset we expect the central diagonal to reach
    // subtract 1 because of the zero index
    // let a_offset: u32 = max![tlen, qlen];
    let a_offset: u32 = tlen;

    // eprintln!("\t a_k {} a_offset {}", a_k, a_offset);

    let max_possible_score = max![
        // longer * mismatch_score
        config.penalties.mismatch as u32 * a_offset,
        // gap_ext * longer + gap_open
        config.penalties.gap_extend as u32 * a_offset + config.penalties.gap_open as u32
    ] as usize;

    let hi: i32 = 0;
    let lo: i32 = 0;

    // Initial conditions
    let wf_set = types::WaveFrontSet {
        i: None,
        d: None,
        m: Some(types::WaveFront::new(hi, lo)),
    };

    let mut all_wavefronts = types::WaveFronts {
        wavefront_set: vec![Some(wf_set)],
        min_k: -(qlen as isize),
        max_k: tlen as isize,
        a_k,
    };

    // score
    let mut score: usize = 0;

    /*
    // Progress bar
    let progress_bar = ProgressBar::new(a_offset as u64);
    let template = "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}]  {pos:>7}/{len:7}  {msg} segments ({eta_precise})";
    let progress_style = ProgressStyle::default_bar()
        .template(template)
        .progress_chars("=> ");
    progress_bar.set_style(progress_style);
    */

    // set score at start to 0
    // unnecessary
    // score ... diagonal
    // all_wavefronts.wavefront_set[0].m.vals[0] = 0;
    if *all_wavefronts
        .get_m_wavefront(score as i32)
        .expect(&format!(
            "[wfa::wf_align] no m-wavefront at score {}",
            score
        ))
        .get_offset(0 as i32)
        .expect(&format!(
            "[wfa::wf_align] no offset on a_k = {a_k} m-wavefront at score {score}"
        ))
        != types::Offset::from([0])
    {
        panic!("[wfa::wf_align] start score should be zero");
    }

    // Print config
    if config.verbosity > 0 {
        eprintln!(
            "Config {{\n\
             {0:two_spaces$}tlen: {1},\n\
             {0:two_spaces$}qlen: {2}\n\
             {0:two_spaces$}a_k: {3}\n\
             {0:two_spaces$}a_offset: {4}\n\
             }}",
            "",
            tlen,
            qlen,
            a_k,
            a_offset,
            two_spaces = 2
        );
    }
    //let mut prev: u64 = 0;
    loop {
        /*
        // update the progress bar
        if let Some(offsets) = all_wavefronts
            .get_m_wavefront(score as i32)
            .and_then(|m_wf| m_wf.get_offset(a_k))
        {
            let m_s_k = offsets.max() as u64;
            let delta = m_s_k - prev;

            progress_bar.inc(delta);
            prev = m_s_k;
        };
         */

        // Extend the current wavefront
        if all_wavefronts.get_m_wavefront(score as i32).is_some() {
            let m_wf_mut: &mut types::WaveFront = &mut all_wavefronts.wavefront_set[score]
                .as_mut()
                .unwrap()
                .m
                .as_mut()
                .unwrap();

            wf_extend(m_wf_mut, match_lambda, &config, score);
        }

        // give up
        if score >= max_possible_score {
            let e = format!("Gave up. {} {}", score, max_possible_score);
            return Err(e);
        }
        // Check whether we have reached the final point
        // Get the m-wavefront with the current score
        if utils::end_reached(
            all_wavefronts.get_m_wavefront(score as i32),
            a_k,
            a_offset,
            config,
        ) {
            break;
        }

        score += 1;

        // TODO: compute the next wavefront
        core::wf_next(&mut all_wavefronts, score, config);
    }

    let cigar = wf_traceback(&all_wavefronts, score, config, traceback_lambda);

    Ok((score, cigar))
}

#[cfg(test)]
mod tests {

    use super::*;

    fn test_config() -> types::Config {
        types::Config {
            adapt: false,
            verbosity: 1,
            penalties: types::Penalties {
                mismatch: 1,
                matches: 0,
                gap_open: 2,
                gap_extend: 1,
            },
        }
    }

    #[test]
    fn test_matches() {
        let query = "ATCGAA".as_bytes();
        let ed_string = "ATC{TA,GA}A";

        let config = test_config();
        let edt = eds::EDT::from_str(ed_string);
        let dt: eds::DT = edt.extract_inelastic();

        let verbosity = config.verbosity;

        let tlen = dt.p();
        let qlen = query.len();

        let mut match_lambda = |v: &mut i32, h: &mut i32, offsets: &mut types::Offset| -> bool {
            if verbosity > 1 {
                eprint!("v ({}, {}) ", h, v);
            }

            if *v < 0 || *h < 0 || *h as usize >= tlen || *v as usize >= qlen {
                return false;
            }

            let text_chars: &Vec<u8> = &dt[*h as usize];
            let query_char: u8 = query[*v as usize];
            let z = text_chars.len();

            let l: usize = offsets.data.len();

            if z > l {
                // copy over
                let prev: i32 = offsets.max();
                *offsets = types::Offset::from_vec(&vec![prev; z]);
            }

            if z < l {
                let furthest: i32 = offsets.max();
                *offsets = types::Offset::from_vec(&vec![furthest; z]);
            }

            //let r = row_chars.iter().map(|c| *c as char).collect::<Vec<char>>();

            if verbosity > 2 {
                // eprint!("offsets {:?} q {} row {:?}", offsets, col_char as char, r);
            }

            let mut found = false;
            let mut increment_once = false;

            for (idx, offset) in offsets.data.iter_mut().enumerate() {
                if text_chars[idx] == query_char {
                    *offset += 1;
                    found = true;

                    if increment_once == false {
                        *v += 1;
                        *h += 1;

                        increment_once = true;
                    }
                }
            }

            if verbosity > 2 {
                // eprint!("found={} {:?}", found, offsets);
                eprintln!();
            }

            found
        };

        let mut traceback_lambda = |_query @ (_q_start, _q_stop): (i32, i32),
                                    _text @ (_t_start, _t_stop): (i32, i32)|
         -> bool { true };

        let (score, cigar) = wf_align(
            tlen as u32,
            qlen as u32,
            &config,
            &mut match_lambda,
            &mut traceback_lambda,
        )
        .expect("[wfa::dwflambda::align::tests::test_matches] Alignment failed");

        assert_eq!("MMMMMM", cigar);
        assert_eq!(0, score);
    }

    // also tests replacement
    #[test]
    fn test_artifact_match() {
        let query = "ATCGAA".as_bytes();
        let ed_string = "ATC{TA,GC}A";

        let config = test_config();
        let edt = eds::EDT::from_str(ed_string);
        let dt: eds::DT = edt.extract_inelastic();

        let verbosity = config.verbosity;

        let tlen = dt.p();
        let qlen = query.len();

        let mut match_lambda = |v: &mut i32, h: &mut i32, offsets: &mut types::Offset| -> bool {
            if verbosity > 4 {
                eprint!("v ({}, {})\n", h, v);
            }

            if *v < 0 || *h < 0 || *h as usize >= tlen || *v as usize >= qlen {
                return false;
            }

            let text_chars: &Vec<u8> = &dt[*h as usize];
            let query_char: u8 = query[*v as usize];
            let z = text_chars.len();

            let l: usize = offsets.offset_count();

            if z > l {
                // copy over
                let prev: i32 = offsets.max();
                *offsets = types::Offset::from_vec(&vec![prev; z]);

                if z > 1 {
                    offsets.abandoned_all_null();
                }
            }

            if z < l {
                let furthest: i32 = offsets.max();
                *offsets = types::Offset::from_vec(&vec![furthest; z]);

                if z == 1 {
                    offsets.abdandoned = None;
                }
            }

            if verbosity > 2 {
                // eprint!("offsets {:?} q {} row {:?}", offsets, col_char as char, r);
            }

            let mut found = false;
            let mut increment_once = false;

            for idx in 0..z {
                if text_chars[idx] == query_char && !offsets.is_abandoned(idx) {
                    offsets.data[idx] += 1;

                    found = true;

                    if increment_once == false {
                        *v += 1;
                        *h += 1;

                        increment_once = true;
                    }
                }

                if text_chars[idx] != query_char {
                    offsets.set_abandon(idx);
                }
            }

            if verbosity > 4 {
                eprint!("\tfound={} {:?}", found, offsets);
                eprintln!();
            }

            found
        };

        let mut traceback_lambda =
            |q @ (q_start, q_stop): (i32, i32), t @ (t_start, t_stop): (i32, i32)| -> bool {
                if q_start < 0 || q_stop < 0 || t_start < 0 || t_stop < 0 {
                    return false;
                }

                let res = (q_start as usize..q_stop as usize)
                    .zip(t_start as usize..t_stop as usize)
                    .fold(true, |acc, (q_index, t_index)| {
                        dt[t_index]
                            .iter()
                            .copied()
                            .any(|t_char| t_char == query[q_index])
                            && acc
                    });

                // eprintln!("{:?} {:?} {}", q, t, res);
                res
            };

        let (score, cigar) = wf_align(
            tlen as u32,
            qlen as u32,
            &config,
            &mut match_lambda,
            &mut traceback_lambda,
        )
        .expect("[wfa::dwflambda::align::tests::test_artifact_match] Alignment failed");

        assert_eq!("MMMMXM", cigar);
        assert_eq!(1, score);
    }

    #[test]
    fn test_snps() {
        let query = "TGGGCACTATCCCTTGTACGTTCGGAGTTTCATATTGTGTATCAAATATATTTATTAG\
                     CTCTTTTGAGCCTGACGAGCTGGGTAG";
        let query = query.as_bytes();
        let ed_string = "TAGGC{TGG,ACT}ATCCCTT{TAA,GTA}{AT,CG}TTCTCA{C,G}TTTC\
                         CA{TGG,ATT}{C,G}TGAATCAAATGTATTTAT{TCGG,TAGG}CT{A,C}TT\
                         TTGAGC{AG,CT}GACTA{GTT,GCT}AGTTAG";

        let config = test_config();
        let edt = eds::EDT::from_str(ed_string);
        let dt: eds::DT = edt.extract_inelastic();

        let verbosity = config.verbosity;

        let tlen = dt.p();
        let qlen = query.len();

        let mut match_lambda = |v: &mut i32, h: &mut i32, offsets: &mut types::Offset| -> bool {
            if verbosity > 4 {
                eprint!("v ({}, {})\n", h, v);
            }

            if *v < 0 || *h < 0 || *h as usize >= tlen || *v as usize >= qlen {
                return false;
            }

            let text_chars: &Vec<u8> = &dt[*h as usize];
            let query_char: u8 = query[*v as usize];
            let z = text_chars.len();

            let l: usize = offsets.offset_count();

            if z > l {
                // copy over
                let prev: i32 = offsets.max();
                *offsets = types::Offset::from_vec(&vec![prev; z]);

                if z > 1 {
                    offsets.abandoned_all_null();
                }
            }

            if z < l {
                let furthest: i32 = offsets.max();
                *offsets = types::Offset::from_vec(&vec![furthest; z]);

                if z == 1 {
                    offsets.abdandoned = None;
                }
            }

            if verbosity > 2 {
                // eprint!("offsets {:?} q {} row {:?}", offsets, col_char as char, r);
            }

            let mut found = false;
            let mut increment_once = false;

            for idx in 0..z {
                if text_chars[idx] == query_char && !offsets.is_abandoned(idx) {
                    offsets.data[idx] += 1;

                    found = true;

                    if increment_once == false {
                        *v += 1;
                        *h += 1;

                        increment_once = true;
                    }
                }

                if text_chars[idx] != query_char {
                    offsets.set_abandon(idx);
                }
            }

            if verbosity > 4 {
                eprint!("\tfound={} {:?}", found, offsets);
                eprintln!();
            }

            found
        };

        let mut traceback_lambda =
            |q @ (q_start, q_stop): (i32, i32), t @ (t_start, t_stop): (i32, i32)| -> bool {
                if q_start < 0 || q_stop < 0 || t_start < 0 || t_stop < 0 {
                    return false;
                }

                let res = (q_start as usize..q_stop as usize)
                    .zip(t_start as usize..t_stop as usize)
                    .fold(true, |acc, (q_index, t_index)| {
                        dt[t_index]
                            .iter()
                            .copied()
                            .any(|t_char| t_char == query[q_index])
                            && acc
                    });

                // eprintln!("{:?} {:?} {}", q, t, res);
                res
            };

        let (score, cigar) = wf_align(
            tlen as u32,
            qlen as u32,
            &config,
            &mut match_lambda,
            &mut traceback_lambda,
        )
        .expect("[wfa::dwflambda::align::tests::test_artifact_match] Alignment failed");

        // assert_eq!("MMMMXM", cigar);
        // assert_eq!(1, score);
    }
}
