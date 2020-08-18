"""
"""

def comp_curr_invest_retro(reg_id, exp, cost_deg):
    """ Compute cost and MDR of investing current year in retrofitting. No
    maintenance.
    
    Parameter:
        reg_id (int): region where retroffiting applied
        exp (Exposures): exposures
        cost_deg (float): in [0,1]: 0 no willingness to pay (cheap), 
            1 maximum willingness to pay (expensive)
    """
    # cost in 0.01 - 0.2 range
    max_cost, min_cost = 0.2, 0.01
    cost_per = (max_cost - min_cost)*cost_deg + min_cost
    print(cost_per)

    # expensive (20%) -> dam_0 damage, cheap (1%) -> dam_f damage
    cost_0, cost_f = 0.2, 0.01
    dam_0, dam_f = 0.7, 0.9
    dam_a = (dam_f - dam_0)/(cost_f-cost_0)
    dam_b = dam_f - dam_a*(cost_f-cost_0)
    dam_per = dam_a*(cost_per-cost_0) + dam_b

    # apply region
    if not reg_id:
        value_area = exp.value.sum()
    else:
        value_area = exp.value[exp.region_id==reg_id].sum()
    
    # percentage of old building in area
    per_old_buildings = 0.35

    # total cost
    return value_area*per_old_buildings*cost_per, dam_per
