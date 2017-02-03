# PhaseField; Small deformations, linear elastic, beam tension (PFL)
AddTest(
    NAME PhaseField_PFL_beam_1e2_tension
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 5e-2 RELTOL 5e-1
    DIFF_DATA
    beam_1e2_out_anal.vtu beam_1e2_out_pcs_0_ts_1500000_t_15.000000.vtu phasefield_anal phasefield
)
# PhaseField; Small deformations, linear elastic, crack propagation (PFL)
AddTest(
    NAME PhaseField_PFL_square_1e3_tension
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_open_line.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-10 RELTOL 1e-10
    DIFF_DATA
    expected_square_c_quad_line_pcs_0_ts_5000_t_0.500000.vtu square_c_quad_line_pcs_0_ts_5000_t_0.500000.vtu phasefield phasefield
)
