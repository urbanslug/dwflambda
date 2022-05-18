use super::types;
use super::utils;

pub fn wf_extend<F>(
    m_wavefront: &mut types::WaveFront,
    match_lambda: &mut F,
    config: &types::Config,
    score: usize,
) where
    F: FnMut(&mut i32, &mut i32, &mut types::Offset) -> bool,
{
    let verbosity = config.verbosity;

    if verbosity > 1 {
        eprintln!("\t[wflambda::wf_extend]");
    }

    // eprintln!("\t\tlo {} hi {}",  m_wavefront.lo, m_wavefront.hi);
    // eprintln!("\t\tscore {}", score);

    for k in m_wavefront.lo..=m_wavefront.hi {
        // let k_index: usize = utils::compute_k_index(m_wavefront.len(), k, m_wavefront.hi);

        // assuming tlen > qlen
        // let m_s_k: i32 = m_wavefront.offsets[k_index];

        // let m_s_k: &types::Offset = m_wavefront.get_offset(k).expect("[wflambda::wf_extend] fail unwrap k={k}");
        // let m_s_k = m_s_k.max();
        // let mut _v = utils::compute_v(m_s_k, k);
        // let mut _h = utils::compute_h(m_s_k, k);

        let offsets: &types::Offset = m_wavefront.foo(k).unwrap();
        let mut vv: i32 = utils::compute_v_new(offsets, k);
        let mut hh: i32 = utils::compute_h_new(offsets, k);

        //let mut v = utils::compute_v(m_s_k, k);
        //let mut h = utils::compute_h(m_s_k, k);

        // eprintln!("\t\t\tvv hh({}, {}) \t\tv h({}, {})", vv, hh, v, h);
        // eprintln!("\t\t\tk {}\toffset {:?}\t\t({}, {})", k, offsets, v, h);

        if verbosity > 5 {
            // eprintln!("\t\t\tk {}\toffset {}\t({}, {})", k, m_s_k, v, h);
            eprintln!("\t {:?} {}", offsets, score);
        }

        // eprintln!("\t\t\tk {}\toffset {}\t({}, {})", k, m_s_k, v, h);
        // vt[v][h] = m_wavefront.vals[k_index] as i32;
        let offsets: &mut types::Offset = m_wavefront.foo(k).unwrap();
        while match_lambda(&mut vv, &mut hh, offsets) {
            if config.verbosity > 6 {
                eprintln!(
                    "\t[wflambda::wf_extend]\n\
                           \t\t({vv} {hh})"
                );
            }
        }
    }
}
