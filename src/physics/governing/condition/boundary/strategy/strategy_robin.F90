submodule(boundary_strategy) strategy_robin
    implicit none
contains

    !> Exchange with a medium, \(q_\mathrm{in}=-h(u)\,(u-u_\infty)\).
    !>
    !> values(1) = h, values(2) = \(u_\infty\). With values(5) = w > 0 the
    !> coefficient switches with the boundary value: values(3) = \(h_b\), the
    !> coefficient below values(4) = \(u_s\), reached over a C1 step of width w,
    !>
    !> \[ h(u)=h_b+(h-h_b)\,S\bigl((u-u_s)/w+\tfrac12\bigr),\qquad S(s)=s^2(3-2s)\ \text{on}\ [0,1], \]
    !>
    !> as a contact whose resistance changes when the surface freezes.
    module subroutine evaluate_robin_bc(self, current_time, u_curr, result)
        implicit none
        class(type_bc_robin), intent(inout) :: self
        real(real64), intent(in) :: current_time
        real(real64), intent(in) :: u_curr
        type(type_bc_result), intent(inout) :: result

        real(real64) :: values(BC_DATA_SIZE)
        real(real64) :: coefficient, slope, fraction

        call self%provider%get_data(current_time, values)
        call result%initialize()

        coefficient = values(1)
        slope = 0.0d0
        if (values(5) > 0.0d0) then
            fraction = (u_curr - values(4))/values(5) + 0.5d0
            if (fraction <= 0.0d0) then
                coefficient = values(3)
            else if (fraction < 1.0d0) then
                coefficient = values(3) + (values(1) - values(3))*fraction**2*(3.0d0 - 2.0d0*fraction)
                slope = (values(1) - values(3))*6.0d0*fraction*(1.0d0 - fraction)/values(5)
            end if
        end if

        result%flux_value = -coefficient*(u_curr - values(2))
        result%flux_derivative = -(coefficient + slope*(u_curr - values(2)))
    end subroutine evaluate_robin_bc

end submodule strategy_robin
