use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::stats_information_theory;

/// Computes the Shannon entropy of a probability distribution.

///

/// Takes a bincode-serialized `Vec<Expr>` representing the probabilities.

/// Returns a bincode-serialized `Expr` representing the entropy.

#[no_mangle]

pub extern "C" fn rssn_bincode_shannon_entropy(
    probs_buf: BincodeBuffer
) -> BincodeBuffer {

    let probs: Option<Vec<Expr>> =
        from_bincode_buffer(&probs_buf);

    if let Some(p) = probs {

        let res = stats_information_theory::shannon_entropy(&p);

        to_bincode_buffer(&res)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Kullback-Leibler divergence (relative entropy) between two probability distributions.

///

/// Takes two bincode-serialized `Vec<Expr>` representing the probability distributions (`p` and `q`).

/// Returns a bincode-serialized `Expr` representing the KL divergence.

#[no_mangle]

pub extern "C" fn rssn_bincode_kl_divergence(
    p_probs_buf: BincodeBuffer,
    q_probs_buf: BincodeBuffer,
) -> BincodeBuffer {

    let p_probs: Option<Vec<Expr>> =
        from_bincode_buffer(
            &p_probs_buf,
        );

    let q_probs: Option<Vec<Expr>> =
        from_bincode_buffer(
            &q_probs_buf,
        );

    if let (Some(p), Some(q)) =
        (p_probs, q_probs)
    {

        match stats_information_theory::kl_divergence(&p, &q) {
            | Ok(res) => to_bincode_buffer(&res),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the cross-entropy between two probability distributions.

///

/// Takes two bincode-serialized `Vec<Expr>` representing the probability distributions (`p` and `q`).

/// Returns a bincode-serialized `Expr` representing the cross-entropy.

#[no_mangle]

pub extern "C" fn rssn_bincode_cross_entropy(
    p_probs_buf: BincodeBuffer,
    q_probs_buf: BincodeBuffer,
) -> BincodeBuffer {

    let p_probs: Option<Vec<Expr>> =
        from_bincode_buffer(
            &p_probs_buf,
        );

    let q_probs: Option<Vec<Expr>> =
        from_bincode_buffer(
            &q_probs_buf,
        );

    if let (Some(p), Some(q)) =
        (p_probs, q_probs)
    {

        match stats_information_theory::cross_entropy(&p, &q) {
            | Ok(res) => to_bincode_buffer(&res),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Gini impurity of a probability distribution.

///

/// Takes a bincode-serialized `Vec<Expr>` representing the probabilities.

/// Returns a bincode-serialized `Expr` representing the Gini impurity.

#[no_mangle]

pub extern "C" fn rssn_bincode_gini_impurity(
    probs_buf: BincodeBuffer
) -> BincodeBuffer {

    let probs: Option<Vec<Expr>> =
        from_bincode_buffer(&probs_buf);

    if let Some(p) = probs {

        let res = stats_information_theory::gini_impurity(&p);

        to_bincode_buffer(&res)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the joint entropy of a joint probability distribution.

///

/// Takes a bincode-serialized `Expr` representing the joint probability distribution.

/// Returns a bincode-serialized `Expr` representing the joint entropy.

#[no_mangle]

pub extern "C" fn rssn_bincode_joint_entropy(
    joint_probs_buf: BincodeBuffer
) -> BincodeBuffer {

    let joint: Option<Expr> =
        from_bincode_buffer(
            &joint_probs_buf,
        );

    if let Some(j) = joint {

        match stats_information_theory::joint_entropy(&j) {
            | Ok(res) => to_bincode_buffer(&res),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the conditional entropy of a joint probability distribution.

///

/// Takes a bincode-serialized `Expr` representing the joint probability distribution.

/// Returns a bincode-serialized `Expr` representing the conditional entropy.

#[no_mangle]

pub extern "C" fn rssn_bincode_conditional_entropy(
    joint_probs_buf: BincodeBuffer
) -> BincodeBuffer {

    let joint: Option<Expr> =
        from_bincode_buffer(
            &joint_probs_buf,
        );

    if let Some(j) = joint {

        match stats_information_theory::conditional_entropy(&j) {
            | Ok(res) => to_bincode_buffer(&res),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the mutual information between two random variables from their joint probability distribution.

///

/// Takes a bincode-serialized `Expr` representing the joint probability distribution.

/// Returns a bincode-serialized `Expr` representing the mutual information.

#[no_mangle]

pub extern "C" fn rssn_bincode_mutual_information(
    joint_probs_buf: BincodeBuffer
) -> BincodeBuffer {

    let joint: Option<Expr> =
        from_bincode_buffer(
            &joint_probs_buf,
        );

    if let Some(j) = joint {

        match stats_information_theory::mutual_information(&j) {
            | Ok(res) => to_bincode_buffer(&res),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}
