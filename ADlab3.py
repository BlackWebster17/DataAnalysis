from spyre import server
import pandas as pd

class SetData(server.App):
    title = "Візуалізація даних NOAA"

    inputs = [
        {
            "type": 'dropdown',
            "label": 'Дані NOAA',
            "options": [{"label": "VCI", "value": "VCI"},
                        {"label": "TCI", "value": "TCI"},
                        {"label": "VHI", "value": "VHI"}],
            "key": 'ticker',
            "action_id": "update_data"
        },
        {
            "type": 'dropdown',
            "label": 'Область',
            "options": [{"label": str(i), "value": str(i)} for i in range(1, 26)],
            "key": 'region',
            "action_id": "update_data"
        },
        {
            "type": 'text',
            "label": 'Рік',
            "key": 'year',
            "value": '',
            "action_id": "update_data"
        },
        {
            "type": 'text',
            "label": 'Початковий тиждень',
            "key": 'start_week',
            "value": '',
            "action_id": "update_data"
        },
        {
            "type": 'text',
            "label": 'Кінцевий тиждень',
            "key": 'end_week',
            "value": '',
            "action_id": "update_data"
        }
    ]

    controls = [{
        "type": "hidden",
        "id": "update_data"
    }]

    tabs = ["Графік", "Таблиця"]

    outputs = [
        {
            "type": "plot",
            "id": "plot",
            "control_id": "update_data",
            "tab": "Графік"
        },
        {
            "type": "table",
            "id": "table_id",
            "control_id": "update_data",
            "tab": "Таблиця",
            "on_page_load": True
        }
    ]

    def getData(self, params):
        data = pd.read_csv("EveryData.csv")
        return data

    def getTable(self, params):
        df = self.getData(params)
        region = int(params['region'])
        year = int(params['year'])
        start_week = int(params['start_week'])
        end_week = int(params['end_week'])
        filtered_df = df[(df['Area'] == region) & (df['Year'] == year) & (df['Week'] >= start_week) & (df['Week'] <= end_week)]
        columns = ['Year', 'Week', params['ticker'], 'Area']
        return filtered_df.loc[:, columns]

    def getPlot(self, params):
        region = int(params['region'])
        year = int(params['year'])
        ticker = params['ticker']
        start_week = int(params['start_week'])
        end_week = int(params['end_week'])
        data = self.getData(params)
        filtered_data = data[(data['Area'] == region) & (data['Year'] == year) & (data['Week'] >= start_week) & (data['Week'] <= end_week)]
        plt_obj = filtered_data.plot(x='Week', y=ticker)
        plt_obj.set_ylabel("Значення")
        plt_obj.set_xlabel("Тиждень")
        plt_obj.set_title(f"Графік для області {region} у {year} для вказаних тижнів")
        plot = plt_obj.get_figure()
        return plot

if __name__ == '__main__':
    app = SetData()
    app.launch()
