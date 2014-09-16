USE lab_proc_db;

DROP TABLE IF EXISTS tag_glossary_tbl;

CREATE TABLE tag_glossary_tbl
(
simple_name varchar(100),
units varchar(100)
);

INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('space_time', 's');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('mass_flow_brush_feeder', 'lb/hr');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('mass_brush_feeder', 'lb');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('Feed_Volts', 'V');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('mass_flow_down_brush', 'L/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('mass_flow_entrainment', 'L/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('mass_flow_argon_tracer', 'L/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('mass_flow_feed_vessel_pressure', 'L/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('mass_flow_methane', 'L/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('steam_flow', 'mL/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('exit_gas_flowrate', 'L/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('flow_rate_exit_gas', 'L/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_biomass_feedrate', 'lb/hr');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_mass_flow_entrainment', 'L/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_mass_argon_tracer', 'L/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_mass_feed_vessel_pressure', 'L/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_mass_methane', 'L/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_steam_HPLC_pump', 'mL/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_primary_quench_HPLC_pump', 'mL/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_secondary_quench_HPLC_pump', 'mL/min');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_steam_waterboiler', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_steam_superheater', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_skin_steam_line_to_lance', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_steam_reactor_entry', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_full_port_body_valve', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_furnace_top_zone_element', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_furnace_middle_zone_element', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_furnace_bottom_zone_element', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_skin_tube_top', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_skin_tube_middle', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_skin_tube_bottom', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_bucket_seal_top', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_bucket_seal_bottom', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_exit_gas', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_ash_knockout', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_heat_tape_ash_knockout', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_process_downstream_ash_KO_horizontal', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_heat_tape_downstream_ash_knockout', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_heat_tape_filter_1', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_heat_tape_between_filters', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_heat_tape_filter_2', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_VLS_1_liquid', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_VLS_2_liquid', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('temp_heat_tape_downstream_second_filter', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('TI_622', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('TI_631', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('TI_632', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('TI_641', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('TI_642', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_temp_skin_tube_top', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_temp_skin_tube_middle', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_temp_skin_tube_bottom', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('HX_330A_SP', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('HX_340A_SP', 'C');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pressure_feeder_vessel', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pressure_entrainment', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pressure_reactor_gas_inlet', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pressure_ash_knockout_vessel', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('output_pressure_ash_knockout', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pressure_product_gas_downstream_filters', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pressure_system_set_by_PCV_to_vent', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pressure_emergency_N2_line', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pressure_water_to_steam_system', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pressure_primary_quench', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pressure_secondary_quench', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pressure_analysis_line', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('setpoint_pressure_ash_knockout', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pp_H2O', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pp_CO2', 'psig');
INSERT INTO tag_glossary_tbl (simple_name, units) VALUES ('pp_Ar', 'psig');









