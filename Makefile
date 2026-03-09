setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py && python analysis/analysis.py

dashboard:
	streamlit run dashboard.py
