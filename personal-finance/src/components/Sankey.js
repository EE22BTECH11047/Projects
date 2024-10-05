import React, { Component } from 'react';
import Highcharts from 'highcharts';
import highchartsSankey from 'highcharts/modules/sankey';
import highchartsExporting from 'highcharts/modules/exporting';
import HighchartsReact from 'highcharts-react-official';

import CONFIG from '../config';

highchartsSankey(Highcharts);
highchartsExporting(Highcharts);

class Sankey extends Component {
  constructor(props) {
    super(props);

    let data = [];

    if (props.items) {
      let assets_sum = 0;
      let liabilities_sum = 0;
      let income_sum = 0;
      let expenses_sum = 0;

      props.items.map(function (item) {
        switch (item.type) {
          case CONFIG.ITEM_TYPE.ASSET:
            const monthly_income = (item.amount * item.apr) / 100 / 12;

            assets_sum += parseInt(monthly_income);

            if (monthly_income > 0) {
              data.push([item.name, 'Income', parseInt(monthly_income)]);
            }
            break;

          case CONFIG.ITEM_TYPE.LIABILITY:
            const monthly_expense = (item.amount * item.apr) / 100 / 12;

            liabilities_sum += parseInt(monthly_expense);

            if (monthly_expense > 0) {
              data.push(['Expenses', item.name, parseInt(monthly_expense)]);
            }
            break;

          case CONFIG.ITEM_TYPE.INCOME:
            income_sum += item.amount;

            if (item.amount > 0) {
              data.push([item.name, 'Income', item.amount]);
            }
            break;

          case CONFIG.ITEM_TYPE.EXPENSE:
            expenses_sum += item.amount;

            if (item.amount > 0) {
              data.push(['Expenses', item.name, item.amount]);
            }
            break;

          default:
        }

        return item;
      });

      income_sum += assets_sum;
      expenses_sum += liabilities_sum;

      if (income_sum < expenses_sum) {
        data.push(['Income', 'Expenses', income_sum]);
        data.push(['Deficit', 'Expenses', expenses_sum - income_sum]);
      }

      if (income_sum >= expenses_sum) {
        data.push([
          'Income',
          'Expenses',
          income_sum - (income_sum - expenses_sum),
        ]);
        data.push(['Income', 'Surplus', income_sum - expenses_sum]);
      }
    }

    this.state = {
      options: {
        chart: {
          backgroundColor: 'transparent',
        },
        title: {
          text: 'Sankey Diagram',
        },
        series: [
          {
            keys: ['from', 'to', 'weight'],
            data: data,
            type: 'sankey',
            name: 'Flow',
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

export default Sankey;
