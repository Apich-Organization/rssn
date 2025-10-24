use vergen_gitcl::{
    BuildBuilder, CargoBuilder, Emitter, GitclBuilder, RustcBuilder, SysinfoBuilder,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut emitter = Emitter::default();

    emitter.add_instructions(&BuildBuilder::all_build()?)?;
    emitter.add_instructions(&CargoBuilder::all_cargo()?)?;
    emitter.add_instructions(&GitclBuilder::all_git()?)?;
    emitter.add_instructions(&RustcBuilder::all_rustc()?)?;
    emitter.add_instructions(&SysinfoBuilder::all_sysinfo()?)?;

    emitter.emit()?;

    Ok(())
}
