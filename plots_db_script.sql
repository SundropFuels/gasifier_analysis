USE lab_proc_db;

DROP TABLE IF EXISTS plots_tbl;

CREATE TABLE plots_tbl
(
plot_id int NOT NULL PRIMARY KEY AUTO_INCREMENT,
type VARCHAR(30),
caption VARCHAR(128),
save_loc VARCHAR(128),
figsize_x int,
figsize_y int,
);

CREATE TABLE timeseries_plots
(
plot_id int NOT NULL PRIMARY KEY,
y_labels


