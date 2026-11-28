#![allow(clippy::wildcard_imports)] // Method declarations share checked object/codec helpers.
use super::*;
use crate::kernels::error_correction::Bch;
use crate::kernels::error_correction::PolyGF256;
use crate::kernels::error_correction::ReedSolomonEncoder;
use crate::kernels::finite_field::PrimeFieldElement;
impl NativeType for Bch {
    const KIND: &'static str = "rssn/native/discrete/bch/v1";
    const MAX_BYTES: usize = 2097152;

    fn validate(&self) -> Result<(), OpError> {
        if self.m() > 16 || self.t() > 32 {
            return Err(error("BCH field/correction budget exceeded"));
        }
        Ok(())
    }
}

impl NativeType for ReedSolomonEncoder {
    const KIND: &'static str = "rssn/native/discrete/reed-solomon-encoder/v1";
    const MAX_BYTES: usize = 65_536;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        std::mem::size_of::<Self>().saturating_add(256_usize.saturating_mul(self.parity_symbols()))
    }
}

object_binding!(
    rs_encoder_new,
    "num_native_rs_encoder_new",
    "crate::kernels::error_correction::ReedSolomonEncoder::new",
    "Prepares a Reed–Solomon encoder for the supplied number of parity bytes. The session-owned encoder retains the generator's feedback table and accepts blocks of at most 255 bytes including parity.",
    true,
    (parity: usize) -> NewObject<ReedSolomonEncoder>,
    |_, parity| Ok(NewObject::normalized(ReedSolomonEncoder::new(parity).map_err(error)?))
);
object_binding!(
    rs_encoder_encode,
    "num_native_rs_encoder_encode",
    "crate::kernels::error_correction::ReedSolomonEncoder::encode",
    "Encodes a tensor of exact bytes using a retained Reed–Solomon encoder. The result begins with the original bytes and ends with the prepared parity count. Input and output have no physical units.",
    true,
    (encoder: Handle<ReedSolomonEncoder>, message: Vec<u8>) -> Vec<u8>,
    |ctx, encoder, message| encoder.with(ctx, |encoder| encoder.encode(&message).map_err(error))
);
object_binding!(
    rs_encoder_parity_symbols,
    "num_native_rs_encoder_parity_symbols",
    "crate::kernels::error_correction::ReedSolomonEncoder::parity_symbols",
    "Returns the exact parity-byte count of a retained Reed–Solomon encoder.",
    true,
    (encoder: Handle<ReedSolomonEncoder>) -> usize,
    |ctx, encoder| encoder.with(ctx, |encoder| Ok(encoder.parity_symbols()))
);
impl NativeType for PolyGF256 {
    const KIND: &'static str = "rssn/native/discrete/gf256-polynomial/v1";
    const MAX_BYTES: usize = 65536;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.0.len(), 4096)?;
        Ok(())
    }
}
impl NativeType for PrimeFieldElement {
    const KIND: &'static str = "rssn/native/discrete/prime-field/v1";
    const MAX_BYTES: usize = 64;

    fn validate(&self) -> Result<(), OpError> {
        if self.modulus != 0 && !crate::kernels::homology::is_prime(self.modulus) {
            return Err(error(
                "finite field modulus must be prime or modulus-agnostic zero",
            ));
        }
        Ok(())
    }
}
object_binding!(bch_new,"num_native_bch_new","crate::kernels::error_correction::Bch::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(m:usize,t:usize)->NewObject<Bch>,|_ctx:&CallCtx<'_>,m:usize,t:usize|{limit(m,16)?;limit(t,32)?;if m<2||t==0{return Err(error("BCH needs m>=2 and t>=1"));}Ok(NewObject::normalized(Bch::new(m,t).map_err(error)?))});
object_binding!(bch_for_data_len,"num_native_bch_for_data_len","crate::kernels::error_correction::Bch::for_data_len","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(bits:usize,t:usize)->NewObject<Bch>,|_ctx:&CallCtx<'_>,bits:usize,t:usize|{limit(bits,65535)?;limit(t,32)?;if t==0{return Err(error("BCH correction capability must be positive"));}Ok(NewObject::normalized(Bch::for_data_len(bits,t).map_err(error)?))});
object_binding!(bch_k,"num_native_bch_k","crate::kernels::error_correction::Bch::k","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Bch>)->usize,|ctx:&CallCtx<'_>,a:Handle<Bch>|{a.with(ctx,|o|Ok(o.k()))});
object_binding!(bch_m,"num_native_bch_m","crate::kernels::error_correction::Bch::m","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Bch>)->usize,|ctx:&CallCtx<'_>,a:Handle<Bch>|{a.with(ctx,|o|Ok(o.m()))});
object_binding!(bch_n,"num_native_bch_n","crate::kernels::error_correction::Bch::n","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Bch>)->usize,|ctx:&CallCtx<'_>,a:Handle<Bch>|{a.with(ctx,|o|Ok(o.n()))});
object_binding!(bch_t,"num_native_bch_t","crate::kernels::error_correction::Bch::t","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Bch>)->usize,|ctx:&CallCtx<'_>,a:Handle<Bch>|{a.with(ctx,|o|Ok(o.t()))});
object_binding!(bch_parity_len,"num_native_bch_parity_len","crate::kernels::error_correction::Bch::parity_len","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Bch>)->usize,|ctx:&CallCtx<'_>,a:Handle<Bch>|{a.with(ctx,|o|Ok(o.parity_len()))});
object_binding!(bch_generator,"num_native_bch_generator","crate::kernels::error_correction::Bch::generator","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Bch>)->Vec<u8>,|ctx:&CallCtx<'_>,a:Handle<Bch>|{a.with(ctx,|o|Ok(o.generator().to_vec()))});
object_binding!(bch_encode,"num_native_bch_encode","crate::kernels::error_correction::Bch::encode","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Bch>,data:Vec<u8>)->Vec<u8>,|ctx:&CallCtx<'_>,a:Handle<Bch>,data:Vec<u8>|{limit(data.len(),65535)?;if data.iter().any(|v|*v>1){return Err(error("BCH data must be binary"));}a.with(ctx,|o|o.encode(&data).map_err(error))});
object_binding!(bch_decode,"num_native_bch_decode","crate::kernels::error_correction::Bch::decode","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Bch>,data:Vec<u8>)->Vec<u8>,|ctx:&CallCtx<'_>,a:Handle<Bch>,data:Vec<u8>|{limit(data.len(),65535)?;if data.iter().any(|v|*v>1){return Err(error("BCH data must be binary"));}a.with(ctx,|o|o.decode(&data).map_err(error))});
object_binding!(bch_correct,"num_native_bch_correct","crate::kernels::error_correction::Bch::correct","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Bch>,data:Vec<u8>)->(Vec<u8>,usize),|ctx:&CallCtx<'_>,a:Handle<Bch>,data:Vec<u8>|{limit(data.len(),65535)?;if data.iter().any(|v|*v>1){return Err(error("BCH data must be binary"));}a.with(ctx,|o|o.correct(&data).map_err(error))});
object_binding!(bch_is_codeword,"num_native_bch_is_codeword","crate::kernels::error_correction::Bch::is_codeword","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Bch>,data:Vec<u8>)->bool,|ctx:&CallCtx<'_>,a:Handle<Bch>,data:Vec<u8>|{limit(data.len(),65535)?;if data.iter().any(|v|*v>1){return Err(error("BCH data must be binary"));}a.with(ctx,|o|Ok(o.is_codeword(&data)))});
object_binding!(poly_gf256_new,"num_native_poly_gf256_new","crate::kernels::error_correction::PolyGF256::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(data:Vec<u8>)->NewObject<PolyGF256>,|_ctx:&CallCtx<'_>,data:Vec<u8>|{limit(data.len(),4096)?;Ok(NewObject::normalized(PolyGF256::new(data)))});
object_binding!(poly_gf256_degree,"num_native_poly_gf256_degree","crate::kernels::error_correction::PolyGF256::degree","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PolyGF256>)->usize,|ctx:&CallCtx<'_>,a:Handle<PolyGF256>|{a.with(ctx,|o|Ok(o.degree()))});
object_binding!(poly_gf256_derivative,"num_native_poly_gf256_derivative","crate::kernels::error_correction::PolyGF256::derivative","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PolyGF256>)->NewObject<PolyGF256>,|ctx:&CallCtx<'_>,a:Handle<PolyGF256>|{a.with(ctx,|o|Ok(NewObject::normalized(o.derivative())))});
object_binding!(poly_gf256_normalize,"num_native_poly_gf256_normalize","crate::kernels::error_correction::PolyGF256::normalize","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PolyGF256>)->NewObject<PolyGF256>,|ctx:&CallCtx<'_>,a:Handle<PolyGF256>|{a.with(ctx,|o|Ok(NewObject::normalized(o.normalize())))});
object_binding!(poly_gf256_eval,"num_native_poly_gf256_eval","crate::kernels::error_correction::PolyGF256::eval","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PolyGF256>,x:u8)->u8,|ctx:&CallCtx<'_>,a:Handle<PolyGF256>,x:u8|{a.with(ctx,|o|Ok(o.eval(x)))});
object_binding!(poly_gf256_scale,"num_native_poly_gf256_scale","crate::kernels::error_correction::PolyGF256::scale","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PolyGF256>,x:u8)->NewObject<PolyGF256>,|ctx:&CallCtx<'_>,a:Handle<PolyGF256>,x:u8|{a.with(ctx,|o|Ok(NewObject::normalized(o.scale(x))))});
object_binding!(poly_gf256_poly_add,"num_native_poly_gf256_poly_add","crate::kernels::error_correction::PolyGF256::poly_add","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PolyGF256>,b:Handle<PolyGF256>)->NewObject<PolyGF256>,|ctx:&CallCtx<'_>,a:Handle<PolyGF256>,b:Handle<PolyGF256>|{a.with(ctx,|o|b.with(ctx,|v|{if o.0.len().checked_mul(v.0.len()).is_none_or(|n|n>1_048_576){return Err(error("GF256 polynomial product work budget exceeded"));}Ok(NewObject::normalized(o.poly_add(v)))}))});
object_binding!(poly_gf256_poly_sub,"num_native_poly_gf256_poly_sub","crate::kernels::error_correction::PolyGF256::poly_sub","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PolyGF256>,b:Handle<PolyGF256>)->NewObject<PolyGF256>,|ctx:&CallCtx<'_>,a:Handle<PolyGF256>,b:Handle<PolyGF256>|{a.with(ctx,|o|b.with(ctx,|v|{if o.0.len().checked_mul(v.0.len()).is_none_or(|n|n>1_048_576){return Err(error("GF256 polynomial product work budget exceeded"));}Ok(NewObject::normalized(o.poly_sub(v)))}))});
object_binding!(poly_gf256_poly_mul,"num_native_poly_gf256_poly_mul","crate::kernels::error_correction::PolyGF256::poly_mul","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PolyGF256>,b:Handle<PolyGF256>)->NewObject<PolyGF256>,|ctx:&CallCtx<'_>,a:Handle<PolyGF256>,b:Handle<PolyGF256>|{a.with(ctx,|o|b.with(ctx,|v|{if o.0.len().checked_mul(v.0.len()).is_none_or(|n|n>1_048_576){return Err(error("GF256 polynomial product work budget exceeded"));}Ok(NewObject::normalized(o.poly_mul(v)))}))});
object_binding!(poly_gf256_poly_div,"num_native_poly_gf256_poly_div","crate::kernels::error_correction::PolyGF256::poly_div","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PolyGF256>,b:Handle<PolyGF256>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<PolyGF256>,b:Handle<PolyGF256>|{a.with(ctx,|o|b.with(ctx,|v|{let(q,r)=o.poly_div(v).map_err(error)?;let q=NativeOutput::encode(NewObject::normalized(q),ctx)?;let r=NativeOutput::encode(NewObject::normalized(r),ctx)?;Ok(Dynamic(api::record(vec![("item0",q),("item1",r)])))}))});
object_binding!(prime_field_new,"num_native_prime_field_new","crate::kernels::finite_field::PrimeFieldElement::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(value:u64,modulus:u64)->NewObject<PrimeFieldElement>,|_ctx:&CallCtx<'_>,value:u64,modulus:u64|{if modulus!=0&&!crate::kernels::homology::is_prime(modulus){return Err(error("finite field modulus must be prime or zero"));}Ok(NewObject::normalized(PrimeFieldElement::new(value,modulus)))});
object_binding!(prime_field_inverse,"num_native_prime_field_inverse","crate::kernels::finite_field::PrimeFieldElement::inverse","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PrimeFieldElement>)->Maybe<NewObject<PrimeFieldElement>>,|ctx:&CallCtx<'_>,a:Handle<PrimeFieldElement>|{a.with(ctx,|o|Ok(Maybe(o.inverse().map(NewObject::normalized))))});
object_binding!(prime_field_is_agnostic,"num_native_prime_field_is_agnostic","crate::kernels::finite_field::PrimeFieldElement::is_agnostic","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PrimeFieldElement>)->bool,|ctx:&CallCtx<'_>,a:Handle<PrimeFieldElement>|{a.with(ctx,|o|Ok(o.is_agnostic()))});
object_binding!(prime_field_pow,"num_native_prime_field_pow","crate::kernels::finite_field::PrimeFieldElement::pow","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PrimeFieldElement>,exponent:u64)->NewObject<PrimeFieldElement>,|ctx:&CallCtx<'_>,a:Handle<PrimeFieldElement>,exponent:u64|{a.with(ctx,|o|{if o.modulus==0{return Err(error("modulus-agnostic power has no finite field"));}Ok(NewObject::normalized(o.pow(exponent)))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &bch_new::BINDING,
            &bch_for_data_len::BINDING,
            &bch_k::BINDING,
            &bch_m::BINDING,
            &bch_n::BINDING,
            &bch_t::BINDING,
            &bch_parity_len::BINDING,
            &bch_generator::BINDING,
            &bch_encode::BINDING,
            &bch_decode::BINDING,
            &bch_correct::BINDING,
            &bch_is_codeword::BINDING,
            &rs_encoder_new::BINDING,
            &rs_encoder_encode::BINDING,
            &rs_encoder_parity_symbols::BINDING,
            &poly_gf256_new::BINDING,
            &poly_gf256_degree::BINDING,
            &poly_gf256_derivative::BINDING,
            &poly_gf256_normalize::BINDING,
            &poly_gf256_eval::BINDING,
            &poly_gf256_scale::BINDING,
            &poly_gf256_poly_add::BINDING,
            &poly_gf256_poly_sub::BINDING,
            &poly_gf256_poly_mul::BINDING,
            &poly_gf256_poly_div::BINDING,
            &prime_field_new::BINDING,
            &prime_field_inverse::BINDING,
            &prime_field_is_agnostic::BINDING,
            &prime_field_pow::BINDING,
        ],
        "numerical_vector",
    )
}
