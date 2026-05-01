#%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

library(geowarp)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 GeoWarp competing model factory ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 model_factory: build a GeoWarp model specification ----
#
# Arguments:
#   df            - data frame with x, y, d columns
#   variable_name - name of the response variable
#   m_formula     - mean model formula (default ~1)
#   mbasis        - use vertical basis functions in mean model
#   vbasis        - use vertical basis functions in variance model
#   vert_warp     - use Bernstein AWU vertical warping (TRUE) or linear (FALSE)
#   geom_warp     - use geometric warping unit
#   warping_mapping - integer vector mapping axes to AWUs
#   covariance    - covariance function name
#   dev_model     - deviation model type: "full", "vert", or "white" (LM)
#   min_depth / max_depth - vertical domain bounds
model_factory <- function(df, variable_name,
                          m_formula = "~1",
                          mbasis = T,
                          vbasis = T,
                          vert_warp = T,
                          geom_warp = T,
                          warping_mapping = c(1L,2L,3L),
                          covariance = c("matern15"),
                          dev_model = "full",
                          min_depth = 0, max_depth = 21){

  ## 2.2 Helper: horizontal domain range ----
  get_horizontal_domains <- function(df) {
    if ('horizontal' %in% colnames(df)) {
      range(df$horizontal)
    } else {
      list(range(df$x), range(df$y))
    }
  }

  ## 2.3 Vertical warping unit ----
  if(vert_warp){
    shared_vertical_warping <- geowarp_bernstein_awu(
      order = (max_depth-min_depth),
      prior = list(type = 'gamma', shape = 1.01, rate = 0.01, lower = 0, upper = Inf)
    )
  }else{
    shared_vertical_warping <- geowarp_linear_awu(
      prior = list(type = 'uniform', shape = 1.01, rate = 0.01, lower = 0, upper = 1000)
    )
  }

  ## 2.4 Geometric warping unit ----
  if(geom_warp){
    geometric_warp <- geowarp_geometric_warping_unit(prior_shape = 4)
  }else{
    geometric_warp <- NULL
  }

  ## 2.5 Deviation process ----
  if(dev_model == "full"){
    deviation_process <- geowarp_deviation_model(
      axial_warping_units = c(
        list(geowarp_linear_awu(
          scaling = 1,
          prior = list(type = 'uniform',shape=0,rate=0,
                       lower = 0,upper = 100))),
        list(geowarp_linear_awu(
          scaling = 1,
          prior = list(type = 'uniform',shape=0,rate=0,
                       lower = 0,upper = 100))),
        list(shared_vertical_warping)
      ),
      axial_warping_unit_mapping = warping_mapping,
      geometric_warping_unit = geometric_warp,
      variance_model =  geowarp_variance_model(
        vertical_basis_functions = vbasis,
        vertical_basis_function_delta = 3,
        vertical_basis_function_boundary_knots = 3L),
      covariance_function = covariance
    )
  }else if(dev_model == "vert"){
    deviation_process <- geowarp_vertical_only_deviation_model(
      axial_warping_unit = shared_vertical_warping,
      variance_model =  geowarp_variance_model(
        vertical_basis_functions = vbasis,
        vertical_basis_function_delta = 5,
        vertical_basis_function_boundary_knots = 3L),
      covariance_function = covariance
    )
  }else{
    deviation_process <- geowarp_white_deviation_model(
      variance_model =  geowarp_variance_model(
        vertical_basis_functions = vbasis,
        vertical_basis_function_delta = 3,
        vertical_basis_function_boundary_knots = 3L)
    )
  }

  ## 2.6 Assemble and return model object ----
  geowarp_model(
    variable = variable_name,
    horizontal_coordinates = c('x', 'y'),
    vertical_coordinate = 'd',
    vertical_domain = c(min_depth, max_depth),
    horizontal_domains = get_horizontal_domains(df),
    nugget_prior = inverse_gamma_quantile_prior(0.001, 1000),
    mean_model = geowarp_mean_model(
      fixed_formula = m_formula,
     vertical_basis_functions = mbasis,
     vertical_basis_function_delta = 3,
    ),
    deviation_model = deviation_process
  )
}
