# Master Thesis Project

### Environment recommendation: 
- Python == 3.9
- R == 4.2 (or earlier version with working Biobase package)

### HOW-TO run the Streamlit App

1. Clone this repository
2. Create new conda environment
3. Install requirements.txt
4. Install Snowphlake
5. Within R, install Biobase, NMF, ggseg3d and htmlwidgets packages 
6. Pass your data by editing *input_data.py* file
7. Run setup file from command line:

>> python app_setup.py

5. Wait 2-3 minutes before Streamlit App is ready. You can track the progress in the cmd/terminal. 

6. Run Streamlit App from command line:

>> streamlit run app.py

