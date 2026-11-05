module module_boundary
    use :: physics_governing_boundary_manager, only: type_bc_manager
    use :: boundary_strategy, only: type_bc_result
    use :: boundary_data_provider, only: BC_DATA_SIZE
    implicit none
    private

    public :: type_bc_manager
    public :: type_bc_result
    public :: BC_DATA_SIZE
end module module_boundary