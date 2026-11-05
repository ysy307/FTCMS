module module_governing
    use :: physics_governing_thermal, only: &
        type_thermal
    use :: physics_governing_hydraulic, only: &
        type_hydraulic
    use :: governing_atmosphere
    use :: module_boundary
    implicit none
    public

end module module_governing
