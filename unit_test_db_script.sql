USE gas_unit_test;

DROP TABLE IF EXISTS gas_data_test_tbl;

CREATE TABLE gas_data_test_tbl
(
timestamp datetime PRIMARY KEY,
counter int,
ME_101 float,
TE_101 float,
TE_102 float,
PT_101 float,
PT_102 float,
PT_103 float,
MFC_101 float,
MFC_102 float,
MFC_103 float,
MFC_104 float,
FE_101 float,
CO_NDIR float,
CO2_NDIR float,
CH4_NDIR float,
H2_GC float,
CO_GC float,
CO2_GC float,
CH4_GC float,
C2H6_GC float,
N2_GC float,
Ar_GC float
);

INSERT INTO gas_data_test_tbl (timestamp, counter, ME_101, TE_101, TE_102, PT_101, PT_102, PT_103, MFC_101, MFC_102, MFC_103, MFC_104, FE_101, CO_NDIR, CO2_NDIR, CH4_NDIR, H2_GC, CO_GC, CO2_GC, CH4_GC, C2H6_GC, N2_GC, Ar_GC) VALUES ('1981-07-06 13:13:12', 0, 3.5, 25.0, 120, 2.51, 1.01, 0.33, 0.0, 7.1, 0.51, 1.01, 51.0, 0.26, 0.09, 0.042, 0.41, 0.26, 0.09, 0.042, 0.01, 0.186, 0.0164);
INSERT INTO gas_data_test_tbl (timestamp, counter, ME_101, TE_101, TE_102, PT_101, PT_102, PT_103, MFC_101, MFC_102, MFC_103, MFC_104, FE_101, CO_NDIR, CO2_NDIR, CH4_NDIR, H2_GC, CO_GC, CO2_GC, CH4_GC, C2H6_GC, N2_GC, Ar_GC) VALUES ('1981-07-06 13:13:13', 1, 3.4, 26.0, 121, 2.5, 1.11, 0.34, 0.0, 7.09, 0.5, 1.0, 51.3, 0.27, 0.089, 0.041, 0.4, 0.27, 0.089, 0.041, 0.01, 0.179, 0.0112);
INSERT INTO gas_data_test_tbl (timestamp, counter, ME_101, TE_101, TE_102, PT_101, PT_102, PT_103, MFC_101, MFC_102, MFC_103, MFC_104, FE_101, CO_NDIR, CO2_NDIR, CH4_NDIR, H2_GC, CO_GC, CO2_GC, CH4_GC, C2H6_GC, N2_GC, Ar_GC) VALUES ('1981-07-06 13:13:14', 2, 3.6, 25.0, 122, 2.55, 1.21, 0.32, 0.0, 7.11, 0.5, 0.97, 50.9, 0.26, 0.088, 0.042, 0.42, 0.26, 0.088, 0.042, 0.011, 0.168, 0.0105);
INSERT INTO gas_data_test_tbl (timestamp, counter, ME_101, TE_101, TE_102, PT_101, PT_102, PT_103, MFC_101, MFC_102, MFC_103, MFC_104, FE_101, CO_NDIR, CO2_NDIR, CH4_NDIR, H2_GC, CO_GC, CO2_GC, CH4_GC, C2H6_GC, N2_GC, Ar_GC) VALUES ('1981-07-06 13:13:15', 3, 3.5, 27.0, 121, 2.46, 1.0, 0.33, 0.0, 7.1, 0.49, 0.99, 50.8, 0.27, 0.089, 0.041, 0.42, 0.27, 0.089, 0.041, 0.01, 0.16, 0.01);
INSERT INTO gas_data_test_tbl (timestamp, counter, ME_101, TE_101, TE_102, PT_101, PT_102, PT_103, MFC_101, MFC_102, MFC_103, MFC_104, FE_101, CO_NDIR, CO2_NDIR, CH4_NDIR, H2_GC, CO_GC, CO2_GC, CH4_GC, C2H6_GC, N2_GC, Ar_GC) VALUES ('1981-07-06 13:13:16', 4, 3.5, 25.0, 122, 2.49, 1.09, 0.34, 0.0, 7.1, 0.5, 1.05, 51.2, 0.28, 0.091, 0.04, 0.42, 0.28, 0.091, 0.04, 0.011, 0.1487, 0.00929);

DROP TABLE IF EXISTS tag_glossary_tbl;

CREATE TABLE tag_glossary_tbl
(
tag varchar(100),
simple_name varchar(100),
units varchar(100)
);

INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('timestamp', 'timestamp', null);
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('ME_101', 'ME_101', 'lb/hr');
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('TE_101', 'TE_101', 'K');
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('TE_102', 'TE_102', 'K');
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('PT_101', 'PT_101', 'lbf/in^2');
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('PT_102', 'PT_102', 'lbf/in^2');
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('PT_103', 'PT_103', 'lbf/in^2');
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('MFC_101', 'MFC_101', 'L/min');
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('MFC_102', 'MFC_102', 'L/min');
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('MFC_103', 'MFC_103', 'L/min');
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('MFC_104', 'MFC_104', 'L/min');
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('FE_101', 'FE_101', 'L/min');
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('CO_NDIR', 'CO_NDIR', null);
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('CO2_NDIR', 'CO2_NDIR', null);
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('CH4_NDIR', 'CH4_NDIR', null);
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('H2_GC', 'H2_GC', null);
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('CO_GC', 'CO_GC', null);
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('CO2_GC', 'CO2_GC', null);
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('CH4_GC', 'CH4_GC', null);
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('C2H6', 'C2H6', null);
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('N2_GC', 'N2_GC', null);
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('Ar_GC', 'Ar_GC', null);
INSERT INTO tag_glossary_tbl (tag, simple_name, units) VALUES ('counter', 'counter', null);
