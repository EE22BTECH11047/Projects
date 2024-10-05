import React, { Component } from 'react';
import Highcharts from 'highcharts';
import highchartsExporting from 'highcharts/modules/exporting';
import HighchartsReact from 'highcharts-react-official';

import CONFIG from '../config';

highchartsExporting(Highcharts);

class Donut extends Component {
  constructor(props) {
    super(props);

    let data = [];

    if (props.items) {
      let expenses_sum = 0;
      let income_sum = 0;

      props.items.map(function (item) {
        switch (item.type) {
          case CONFIG.ITEM_TYPE.ASSET:
            const monthly_income = (item.amount * item.apr) / 100 / 12;

            income_sum += monthly_income;
            break;

          case CONFIG.ITEM_TYPE.LIABILITY:
            const monthly_expense = (item.amount * item.apr) / 100 / 12;

            expenses_sum += monthly_expense;
            break;

          case CONFIG.ITEM_TYPE.INCOME:
            income_sum += item.amount;
            break;

          case CONFIG.ITEM_TYPE.EXPENSE:
            expenses_sum += item.amount;
            break;

          default:
        }

        return item;
      });

      data.push(['Income', income_sum]);
      data.push(['Expenses', expenses_sum]);
    }

    this.state = {
      options: {
        chart: {
          plotBackgroundColor: null,
          plotBorderWidth: 0,
          plotShadow: false,
          backgroundColor: 'transparent',
        },
        title: {
          text: 'Health',
          align: 'center',
          verticalAlign: 'middle',
          y: 60,
        },
        tooltip: {
          pointFormat: '{series.name}: <b>{point.percentage:.1f}%</b>',
        },
        plotOptions: {
          pie: {
            colors: ['#4582EC', '#d9534f'],
            dataLabels: {
              enabled: true,
              distance: -50,
              style: {
                fontWeight: 'bold',
                color: 'white',
              },
            },
            startAngle: -90,
            endAngle: 90,
            center: ['50%', '75%'],
            size: '110%',
          },
        },
        series: [
          {
            type: 'pie',
            name: 'Balance',
            innerSize: '50%',
            data: data,
          },
        ],
      },
    };
  }

  render() {
    return (
      <HighchartsReact highcharts={Highcharts} options={this.state.options} />
    );
  }
}

export default Donut;
