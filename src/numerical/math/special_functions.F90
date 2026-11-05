!> Numerical special functions.
module numerical_special_functions
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none
    private

    integer(int32), parameter :: beta_max_iterations = 200

    public :: type_regularized_incomplete_beta

    !> Regularized incomplete-beta evaluator for fixed \(a,b\).
    type :: type_regularized_incomplete_beta
        private
        real(real64) :: a = 0.0d0
        real(real64) :: b = 0.0d0
        real(real64) :: log_normalization = 0.0d0
        real(real64) :: direct_even_factor(beta_max_iterations) = 0.0d0
        real(real64) :: direct_odd_factor(beta_max_iterations) = 0.0d0
        real(real64) :: swapped_even_factor(beta_max_iterations) = 0.0d0
        real(real64) :: swapped_odd_factor(beta_max_iterations) = 0.0d0
        logical :: initialized = .false.
    contains
        procedure, pass(self), public :: initialize
        procedure, pass(self), public :: evaluate
        procedure, nopass, private :: beta_continued_fraction
        procedure, pass(self), public :: get_log_normalization
    end type type_regularized_incomplete_beta

contains

    !> Cache the \(\ln\Gamma\) normalization for fixed \(a,b\).
    subroutine initialize(self, a, b)
        implicit none
        class(type_regularized_incomplete_beta), intent(inout) :: self
        real(real64), intent(in) :: a
        real(real64), intent(in) :: b

        integer(int32) :: iteration, twice_iteration
        real(real64) :: iteration_real, twice_iteration_real

        self%a = a
        self%b = b
        self%log_normalization = 0.0d0
        self%initialized = a > 0.0d0 .and. b > 0.0d0
        if (.not. self%initialized) return

        self%log_normalization = log_gamma(a + b) - log_gamma(a) - log_gamma(b)

        do iteration = 1, beta_max_iterations
            twice_iteration = 2 * iteration
            iteration_real = real(iteration, real64)
            twice_iteration_real = real(twice_iteration, real64)

            self%direct_even_factor(iteration) = iteration_real * (b - iteration_real) / &
                ((a - 1.0d0 + twice_iteration_real) * (a + twice_iteration_real))
            self%direct_odd_factor(iteration) = -(a + iteration_real) * (a + b + iteration_real) / &
                ((a + twice_iteration_real) * (a + 1.0d0 + twice_iteration_real))

            self%swapped_even_factor(iteration) = iteration_real * (a - iteration_real) / &
                ((b - 1.0d0 + twice_iteration_real) * (b + twice_iteration_real))
            self%swapped_odd_factor(iteration) = -(b + iteration_real) * (a + b + iteration_real) / &
                ((b + twice_iteration_real) * (b + 1.0d0 + twice_iteration_real))
        end do
    end subroutine initialize

    !> Compute the regularized incomplete beta function \(I_x(a,b)\).
    pure subroutine evaluate(self, x, result, converged)
        implicit none
        class(type_regularized_incomplete_beta), intent(in) :: self
        real(real64), intent(in) :: x
        real(real64), intent(inout) :: result
        logical, intent(inout), optional :: converged

        real(real64) :: beta_factor
        real(real64) :: continued_fraction
        logical :: fraction_converged

        if (present(converged)) converged = .true.
        if (x <= 0.0d0) then
            result = 0.0d0
            return
        end if
        if (x >= 1.0d0) then
            result = 1.0d0
            return
        end if
        if (.not. self%initialized) then
            result = 0.0d0
            if (present(converged)) converged = .false.
            return
        end if

        beta_factor = exp(self%log_normalization + self%a * log(x) + self%b * log(1.0d0 - x))
        if (x < (self%a + 1.0d0) / (self%a + self%b + 2.0d0)) then
            call self%beta_continued_fraction(self%a + self%b, self%a + 1.0d0, &
                                               self%direct_even_factor, self%direct_odd_factor, &
                                               x, continued_fraction, fraction_converged)
            result = beta_factor * continued_fraction / self%a
        else
            call self%beta_continued_fraction(self%a + self%b, self%b + 1.0d0, &
                                               self%swapped_even_factor, self%swapped_odd_factor, &
                                               1.0d0 - x, continued_fraction, fraction_converged)
            result = 1.0d0 - beta_factor * continued_fraction / self%b
        end if
        result = min(1.0d0, max(0.0d0, result))
        if (present(converged)) converged = fraction_converged
    end subroutine evaluate

    !> Evaluate the incomplete-beta continued fraction using modified Lentz iteration.
    pure subroutine beta_continued_fraction(a_plus_b, a_plus_one, even_factor, odd_factor, x, result, converged)
        implicit none
        real(real64), intent(in) :: a_plus_b
        real(real64), intent(in) :: a_plus_one
        real(real64), intent(in) :: even_factor(beta_max_iterations)
        real(real64), intent(in) :: odd_factor(beta_max_iterations)
        real(real64), intent(in) :: x
        real(real64), intent(inout) :: result
        logical, intent(inout) :: converged

        real(real64), parameter :: convergence_tolerance = 3.0d-15
        real(real64), parameter :: minimum_fraction = 1.0d-300
        real(real64) :: coefficient, numerator, denominator, delta
        integer(int32) :: iteration

        converged = .false.
        numerator = 1.0d0
        denominator = 1.0d0 - a_plus_b * x / a_plus_one
        if (abs(denominator) < minimum_fraction) denominator = minimum_fraction
        denominator = 1.0d0 / denominator
        result = denominator

        do iteration = 1, beta_max_iterations
            coefficient = even_factor(iteration) * x
            denominator = 1.0d0 + coefficient * denominator
            if (abs(denominator) < minimum_fraction) denominator = minimum_fraction
            numerator = 1.0d0 + coefficient / numerator
            if (abs(numerator) < minimum_fraction) numerator = minimum_fraction
            denominator = 1.0d0 / denominator
            result = result * denominator * numerator

            coefficient = odd_factor(iteration) * x
            denominator = 1.0d0 + coefficient * denominator
            if (abs(denominator) < minimum_fraction) denominator = minimum_fraction
            numerator = 1.0d0 + coefficient / numerator
            if (abs(numerator) < minimum_fraction) numerator = minimum_fraction
            denominator = 1.0d0 / denominator
            delta = denominator * numerator
            result = result * delta
            if (abs(delta - 1.0d0) < convergence_tolerance) then
                converged = .true.
                return
            end if
        end do
    end subroutine beta_continued_fraction

    !> \(-\ln B(a,b)\), the normalization of \(I_x(a,b)\).
    !>
    !> \[ \frac{\mathrm{d}I_x}{\mathrm{d}x}=\exp\bigl(-\ln B(a,b)\bigr)\,x^{a-1}(1-x)^{b-1} \]
    !>
    !> Assumptions: `initialize` was called; zero otherwise.
    !> Numerical guarantees: the value `initialize` cached from \(\ln\Gamma\).
    !> Computational complexity: \(O(1)\).
    !> Failure behavior: none.
    pure function get_log_normalization(self) result(log_normalization)
        implicit none
        class(type_regularized_incomplete_beta), intent(in) :: self
        !> \(-\ln B(a,b)\)
        real(real64) :: log_normalization

        log_normalization = self%log_normalization
    end function get_log_normalization

end module numerical_special_functions
