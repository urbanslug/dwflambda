use super::types;
use super::utils::{self, backtrace as backtrace_utils};

pub fn wf_traceback<G>(
    all_wavefronts: &types::WaveFronts,
    score: usize,
    config: &types::Config,
    traceback_lambda: &mut G,
) -> String
where
    G: FnMut((i32, i32), (i32, i32)) -> bool,
{
    if config.verbosity > 0 {
        eprintln!("\n\t[wfa::wf_backtrace]");
    }

    let x: i32 = config.penalties.mismatch;
    let o: i32 = config.penalties.gap_open;
    let e: i32 = config.penalties.gap_extend;

    let mut cigar = String::new();

    // start with central diagonal
    let mut k = all_wavefronts.a_k as i32;

    // start at the furthest offset on the m-wavefront i.e. the end of the alignment
    let m_wf = all_wavefronts.get_m_wavefront(score as i32).unwrap();
    // let wave_length = m_wf.len();
    // let hi = m_wf.hi;
    // let k_index = utils::compute_k_index(wave_length, k, hi);

    // offset
    let m_s_k: &types::Offset = m_wf
        .get_offset(k)
        .expect("[wflambda::wf_extend] fail unwrap k={k}");
    // let m_s_k: types::Offset = m_wf.offsets[k_index];
    let offsets: &types::Offset = m_wf.bar(k).unwrap();

    let mut v = utils::compute_v_new(offsets, k);
    let mut h = utils::compute_h_new(offsets, k);

    if config.verbosity > 5 {
        eprintln!("\t\t({}, {})", v, h);
    }

    let mut backtrace_op = types::BacktraceOperation::MatchMismatch;

    let mut s = score as i32;
    let mut offset: types::Offset = m_s_k.clone();
    let mut offset: i32 = m_s_k.max();

    // eprintln!("score {}", score);

    while v > 0 && h > 0 && s > 0 {
        // compute scores
        let gap_open_score: i32 = s - o - e;
        let gap_extend_score: i32 = s - e;
        let mismatch_score: i32 = s - x;

        if config.verbosity > 4 {
            eprintln!(
                "\t\tscore: {} \n\
                       \t\tOperation: {:?} \n\
                       \t\t{{\n\
                       \t\t\tg_o: {} \n\
                       \t\t\tg_e: {} \n\
                       \t\t\tx: {} \n\
                       \t\t}}\
                       ",
                s, backtrace_op, gap_open_score, gap_extend_score, mismatch_score
            );
        }

        let del_ext: Option<i32> = if backtrace_op == types::BacktraceOperation::Insertion {
            None
        } else {
            backtrace_utils::backtrace_deletion_extend_offset(all_wavefronts, gap_extend_score, k)
        };

        let del_open: Option<i32> = if backtrace_op == types::BacktraceOperation::Insertion {
            None
        } else {
            backtrace_utils::backtrace_deletion_open_offset(all_wavefronts, gap_open_score, k)
        };

        let ins_ext: Option<i32> = if backtrace_op == types::BacktraceOperation::Deletion {
            None
        } else {
            backtrace_utils::backtrace_insertion_extend_offset(all_wavefronts, gap_extend_score, k)
        };

        let ins_open: Option<i32> = if backtrace_op == types::BacktraceOperation::Deletion {
            None
        } else {
            backtrace_utils::backtrace_insertion_open_offset(all_wavefronts, gap_open_score, k)
        };

        let misms: Option<i32> = if backtrace_op != types::BacktraceOperation::MatchMismatch {
            None
        } else {
            backtrace_utils::backtrace_mismatch_offset(all_wavefronts, mismatch_score, k)
        };

        // Compute maximum offset
        let max_all: Option<i32> = vec![del_ext, del_open, ins_ext, ins_open, misms]
            .into_iter()
            .max()
            .unwrap();

        if config.verbosity > 4 {
            let res = vec![del_ext, del_open, ins_ext, ins_open, misms];
            eprintln!(
                "\t\tdel_ext, del_open, ins_ext, ins_open, misms\n\
                       \t\tops {:?} \n\
                       \t\toffset {:?} \n\
                       \t\tmax_all {:?} \n\
                       \t\tbacktrace_op {:?}",
                res, offset, max_all, backtrace_op
            );
        }

        // Traceback Matches
        if max_all.is_some()
            && backtrace_op == types::BacktraceOperation::MatchMismatch
            && offset >= max_all.unwrap()
        {
            let num_matches = (offset - max_all.unwrap()) as u32;
            backtrace_utils::wflambda_backtrace_matches_check(
                &mut offset,
                &mut cigar,
                num_matches,
                k,
                traceback_lambda,
            );

            offset = max_all.unwrap();
        }

        if max_all == del_ext {
            // Extend a deletion
            cigar.push('D');
            // Update state
            s = gap_extend_score;
            k += 1;
            backtrace_op = types::BacktraceOperation::Deletion;
        } else if max_all == del_open {
            // Open a deletion
            cigar.push('D');
            // Update state
            s = gap_open_score;
            k += 1;
            backtrace_op = types::BacktraceOperation::MatchMismatch;
        } else if max_all == ins_ext {
            // Extend an insertion
            cigar.push('I');
            // Update state
            s = gap_extend_score;
            k -= 1;
            // offset.data.iter_mut().for_each(|offset| *offset -= 1);
            offset -= 1;
            backtrace_op = types::BacktraceOperation::Insertion;
        } else if max_all == ins_open {
            // Add Insertion
            cigar.push('I');
            // Update state
            s = gap_open_score;
            k -= 1;
            // offset.data.iter_mut().for_each(|offset| *offset -= 1);
            offset -= 1;
            backtrace_op = types::BacktraceOperation::MatchMismatch;
        } else if max_all == misms {
            // Add Mismatch
            cigar.push('X');

            // Update state
            s = mismatch_score;
            // offset.data.iter_mut().for_each(|offset| *offset -= 1);
            offset -= 1;
        } else {
            panic!("Backtrace error: No link found during backtrace");
        }

        v = utils::compute_v(offset, k);
        h = utils::compute_h(offset, k);

        if config.verbosity > 5 {
            eprintln!("\t\t({}, {}) s {}", v, h, s);
        }
    }

    // reached the end of one or both of the sequences
    if s == 0 {
        // backtrace matches check
        let num_matches = offset as u32;
        backtrace_utils::wflambda_backtrace_matches_check(
            &mut offset,
            &mut cigar,
            num_matches,
            k,
            traceback_lambda,
        );
    } else {
        // add indels
        while v > 0 {
            cigar.push('D');
            v -= 1;
        }

        while h > 0 {
            cigar.push('I');
            h -= 1;
        }
    }

    let reversed_cigar = cigar.chars().rev().collect::<String>();
    reversed_cigar
}
